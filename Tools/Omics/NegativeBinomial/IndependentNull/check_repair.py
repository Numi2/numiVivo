#!/usr/bin/env python3
"""Verify the follow-up changes only diagnostic-dependent test availability and BH."""
import argparse, gzip, hashlib, json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--available',action='store_true')
a=p.parse_args();root=a.root;results=[]
load=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
for number in range(1,10):
    for method in ['wald','lrt']:
        name=f'{number:02d}-{method}';old=root/'native'/name;new=root/'native-fixed'/name
        if not (new/'receipt.json').exists():assert a.available;continue
        for directory in [old,new]:
            assert hashlib.sha256((directory/'output.json.gz').read_bytes()).hexdigest()==load(directory/'receipt.json')['outputSHA256']
        before=load(old/'output.json.gz')['result'];after=load(new/'output.json.gz')['result']
        assert before['design']==after['design'] and before['request']==after['request'] and before['evidence']==after['evidence']
        assert before['negativeBinomial']['trend']==after['negativeBinomial']['trend']
        old_d=before['negativeBinomial']['features'];new_d=after['negativeBinomial']['features'];expanded=[]
        for i,(b,c) in enumerate(zip(before['features'],after['features'])):
            assert old_d[i].get('finalFit')==new_d[i].get('finalFit')
            assert old_d[i].get('geneWiseDispersion')==new_d[i].get('geneWiseDispersion')
            assert old_d[i].get('finalDispersion')==new_d[i].get('finalDispersion')
            if b['status']==c['status']:
                assert {k:v for k,v in b.items() if k!='adjustedPValue'}=={k:v for k,v in c.items() if k!='adjustedPValue'}
            else:
                assert b['status']=='numericalFailure'
                assert 'influence unavailable on unit-leverage design' in old_d[i].get('error','')
                assert old_d[i]['finalFit'].get('cooksDistances') is None
                assert c['status'] in ['tested','dispersionBoundary']
                expanded.append(c['featureID'])
        results.append(dict(case=name,originalTested=before['testedFeatures'],repairedTested=after['testedFeatures'],
            availabilityChangedFeatures=expanded,identicalCountsDesignNormalizationTrendAndFits=True,existingRawProbabilitiesExactlyEqual=True))
summary=dict(status='passed' if len(results)==18 else 'partial',cases=len(results),results=results)
(root/('repair-check-available.json' if a.available else 'repair-check.json')).write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
print(summary['status'],'cases',len(results),'newly tested',sum(r['repairedTested']-r['originalTested'] for r in results))
