#!/usr/bin/env python3
"""Summarize the complete frozen family, retaining failures and coverage."""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
protocol=json.loads((root/'protocol.json').read_text());native=json.loads((root/'native-complete.json').read_text());reference=json.loads((root/'reference-complete.json').read_text())
assert len(native['runs'])==len(protocol['cases'])==50 and all(r['exitCode']==0 for r in native['runs'])
assert {r['case'] for r in reference}=={r['id'] for r in protocol['cases'] if r['family']=='treatment'}
checks=[];null=[];treatment=[]
for rec in protocol['cases']:
 c=json.loads((root/rec['id']/'check.json').read_text());assert c['case']==rec['id']
 assert c['checkerSHA256']==hashlib.sha256(Path(__file__).with_name('check.py').read_bytes()).hexdigest()
 run=next(r for r in native['runs'] if r['case']==rec['id']);assert c['outputSHA256']==run['outputSHA256'] and c['sourceSHA256']==rec['sourceSHA256']
 checks.append(c)
 if rec['family']=='treatment':
  ref=next(r for r in reference if r['case']==rec['id'])
  treatment.append(dict(case=rec['id'],originalTested=c['oldTestedGenes'],lrtTested=c['testedGenes'],waldBHCalls=c['originalBHCalls'],lrtBHCalls=c['likelihoodRatioBHCalls'],numericalStatus=c['status'],referenceStatus=ref.get('status','failed'),referenceErrors=len(ref.get('errors',[]))))
for study in ['kang','hagai']:
 for policy in ['default','active']:
  cs=[c for c in checks if c['case'].startswith(study+'-null-') and c['case'].endswith('-'+policy)];assert len(cs)==10
  null.append(dict(study=study,policy=policy,testedRange=[min(c['testedGenes'] for c in cs),max(c['testedGenes'] for c in cs)],waldBHCalls=sum(c['originalBHCalls'] for c in cs),lrtBHCalls=sum(c['likelihoodRatioBHCalls'] for c in cs),lrtSplitsWithCalls=sum(c['likelihoodRatioBHCalls']>0 for c in cs),unavailable=sum(len(c['unavailable']) for c in cs)))
failures=[dict(case=c['case'],failures=c['failures'],unavailable=c['unavailable']) for c in checks if c['status']!='passed']
result=dict(status='completed' if not failures and all(r.get('status')=='passed' for r in reference) else 'completed-with-retained-failures',cases=len(checks),testedGeneAnalyses=sum(c['testedGenes'] for c in checks),null=null,treatment=treatment,numericalFailures=failures,referenceFailures=[r for r in reference if r.get('status')!='passed'],maximumErrors={k:max(c['metrics'][k] for c in checks) for k in checks[0]['metrics']},referenceMaximumStatisticDifference=max(r.get('maximumStatisticDifference',0) for r in reference),referenceMaximumPValueDifference=max(r.get('maximumPValueDifference',0) for r in reference),referencePrecisionRefinements=sum(len(c.get('precisionRefinements',[])) for c in checks),protocolSHA256=protocol['protocolSHA256'],binarySHA256=native['binarySHA256'],qualification='Previously inspected data; overlapping sham splits; fixed-dispersion asymptotic LRT; no independent experiment, dispersion-uncertainty adjustment, power or general FDR-control claim')
if (root/'reference-tight-complete.json').exists():
 tight=json.loads((root/'reference-tight-complete.json').read_text());assert {r['case'] for r in tight}=={r['case'] for r in reference}
 changes=[]
 def calls(rows):
  ps=np.asarray([r['pValue'] for r in rows]);order=np.argsort(ps);q=np.empty(len(ps));q[order]=np.minimum(1,np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1]);return {rows[i]['featureIndex'] for i in range(len(rows)) if q[i]<.05}
 for r in tight:
  d=root/r['case'];old=json.loads(gzip.decompress((d/'r-output.json.gz').read_bytes()))['results'];new=json.loads(gzip.decompress((d/'r-output-tight.json.gz').read_bytes()))['results']
  assert {x['featureIndex'] for x in old}=={x['featureIndex'] for x in new}
  changes.append(dict(case=r['case'],changedBHCallIndices=sorted(calls(old)^calls(new))))
 result['tightReference']=dict(cases=len(tight),passed=sum(r.get('status')=='passed' for r in tight),failures=[r for r in tight if r.get('status')!='passed'],maximumStatisticDifference=max(r.get('maximumStatisticDifference',0) for r in tight),maximumPValueDifference=max(r.get('maximumPValueDifference',0) for r in tight),BHChanges=changes,qualification='Supplemental stopping sensitivity; original default-reference failures retained')
if (root/'product/complete.json').exists():
 product=json.loads((root/'product/complete.json').read_text());lrt=json.loads(gzip.decompress((root/'product/likelihoodRatio/report.json.gz').read_bytes()))['contrasts'][0]
 measurement=json.loads(gzip.decompress((root/'kang-default/native.json.gz').read_bytes()))
 assert lrt['features']==[g['feature'] for g in measurement['genes']]
 assert [d.get('likelihoodRatioFit') for d in lrt['negativeBinomial']['features']]==[g.get('likelihoodRatio') for g in measurement['genes']]
 result['product']=dict(**product,harnessFeaturesAndNullFitsExactlyEqual=True)
(root/'summary.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,indent=2))
