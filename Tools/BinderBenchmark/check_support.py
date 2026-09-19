#!/usr/bin/env python3
"""Check actual policy-gated native bundles against already executed legacy folds.
The caller supplies an explicit policy; this never searches thresholds for a better score.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def load(path): return json.loads(Path(path).read_text())


def check(campaign: Path, binary: Path, policy: Path, output: Path):
    output.mkdir(exist_ok=False)
    pol = load(policy)
    summary = []
    for plan_path in sorted((campaign/'plans').glob('*.json')):
        stem=plan_path.stem
        assay=stem.split('-',1)[0]
        imported=campaign/f'input-{assay}'
        data=load(imported/'imported.json')['dataset']
        plan=load(plan_path)
        out=output/stem
        result=subprocess.run([str(binary),'binder-evaluate-supported',str(imported),str(plan_path),str(policy),str(out)],capture_output=True,text=True)
        if result.returncode not in (0,2): raise AssertionError((result.returncode,result.stderr))
        subprocess.run([str(binary),'binder-verify',str(out)],check=True,capture_output=True)
        report=load(out/'report.json'); support=report['support']
        test_groups={r['leakageGroup'] for r in data['records'] if r['target'] in plan['testTargets']}
        required=set(plan['modelFeatures']+[plan['baselineFeature']])
        train=[r for r in data['records'] if r['target'] in plan['trainingTargets'] and r['leakageGroup'] not in test_groups
               and r['outcome'] in ('binder','nonBinder') and required<=set(r['features'])]
        assert support['trainingIDs']==sorted(r['id'] for r in train)
        def counts(rows):
            groups={g:[r for r in rows if r['leakageGroup']==g] for g in {r['leakageGroup'] for r in rows}}
            classes=[{r['outcome'] for r in rs} for rs in groups.values()]
            pos=sum(r['outcome']=='binder' for r in rows)
            return {'candidates':len(rows),'positiveRows':pos,'negativeRows':len(rows)-pos,
                    'positiveGroups':sum(c=={'binder'} for c in classes),
                    'negativeGroups':sum(c=={'nonBinder'} for c in classes),'mixedGroups':sum(len(c)==2 for c in classes)}
        total=counts(train)
        by={t:counts([r for r in train if r['target']==t]) for t in {r['target'] for r in train}}
        assert support['overall']==total and support['byTarget']==by
        eligible=(len(train)>=4 and total['positiveGroups']>=pol['minimumPositiveGroups']
                  and total['negativeGroups']>=pol['minimumNegativeGroups'] and len(by)>=pol['minimumTrainingTargets']
                  and sum(c['positiveGroups']>0 and c['negativeGroups']>0 for c in by.values())>=pol['minimumTargetsWithBothClasses']
                  and len(plan['modelFeatures'])<=pol['maximumFeatures'])
        assert support['eligibleForExperimentalFit']==eligible
        assert result.returncode==(0 if eligible else 2)
        if eligible:
            assert report['fittedReport']==load(campaign/f'result-{stem}'/'report.json')
        else:
            assert 'fittedReport' not in report and report['disposition']=='retainFixedRanking'
        summary.append({'fold':stem,'eligibleForExperimentalFit':eligible,'counts':total,'unmetRequirements':support['unmetRequirements']})
    if len(summary)!=6: raise AssertionError('expected six existing campaign folds')
    data={'schemaVersion':1,'scope':'explicit example policy on reused development data; not a calibrated support rule',
          'policy':pol,'sourceSHA256':load(campaign/'summary.json')['sourceSHA256'],
          'implementationSHA256':hashlib.sha256(binary.read_bytes()).hexdigest(),
          'independentCountsAndRetainedLegacyFits':'passed','results':summary}
    (output/'summary.json').write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print(json.dumps(data,sort_keys=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['campaign','binary','policy','output']: p.add_argument(name,type=Path)
    a=p.parse_args();check(a.campaign.resolve(),a.binary.resolve(),a.policy.resolve(),a.output.absolute())
