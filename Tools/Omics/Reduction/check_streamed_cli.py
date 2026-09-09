#!/usr/bin/env python3
"""Full-cohort streamed PCA publication, reconstruction and bounded failures."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--prepared',type=Path,required=True)
p.add_argument('--mapping',type=Path)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
source=a.prepared/'prepared.h5ad'
mapping=json.loads((a.mapping or a.prepared/'mapping.json').read_text())
plan=dict(schemaVersion=1,mapping=mapping,contrasts=[],reduction={})
plan_path=a.out/'plan.json';plan_path.write_text(json.dumps(plan,indent=2)+'\n')
commands=[]
def run(arguments,expected=None):
    start=time.monotonic()
    result=subprocess.run(['/usr/bin/time','-l',str(a.binary.resolve()),*map(str,arguments)],capture_output=True,text=True,timeout=600)
    i=len(commands);(a.out/f'{i:02}.stdout').write_text(result.stdout);(a.out/f'{i:02}.stderr').write_text(result.stderr)
    def measurement(label):
        return next((int(line.split()[0]) for line in result.stderr.splitlines() if label in line),None)
    commands.append(dict(arguments=list(map(str,arguments)),exitCode=result.returncode,seconds=time.monotonic()-start,
        maximumResidentBytes=measurement('maximum resident set size'),peakFootprintBytes=measurement('peak memory footprint'),expectedFailure=expected))
    (a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    assert result.returncode==0 if expected is None else result.returncode>0 and expected in result.stderr.lower(),result.stderr
    return result
bundle=a.out/'bundle'
run(['singlecell-h5ad-pseudobulk',source,'--plan',plan_path,'--output',bundle])
run(['singlecell-h5ad-pseudobulk-verify',bundle])
repeat=a.out/'repeat'
run(['singlecell-h5ad-pseudobulk',source,'--plan',plan_path,'--output',repeat])
assert (bundle/'receipt.json').read_bytes()==(repeat/'receipt.json').read_bytes()
assert sorted(p.name for p in bundle.iterdir())==['original.h5ad','plan.json','receipt.json','report.json']
for name,options,error in [('cache',{'maximumCacheBytes':16},'cache-byte budget'),
                           ('work',{'maximumEntryVisits':1},'entry-visit budget'),
                           ('basis',{'pca':{'maximumBasis':20}},'residual')]:
    invalid=dict(plan,reduction=options);path=a.out/f'insufficient-{name}.json'
    path.write_text(json.dumps(invalid,indent=2)+'\n');destination=a.out/f'rejected-{name}'
    run(['singlecell-h5ad-pseudobulk',source,'--plan',path,'--output',destination],error)
    assert not destination.exists()
assert not list(a.out.glob('.numivivo-stream-*'))
summary=dict(status='passed',binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
    checks=['full-source streamed PCA publication','native reconstruction','exact repeat',
            'scratch cache absent from published bundle','cache budget rejection','entry visit budget rejection',
            'unconverged PCA rejection','staging cleanup after failure'],commands=commands,
    qualification='CPU workflow and storage/numerical checks; not biological integration or million-cell performance')
(a.out/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='commands'},indent=2))
