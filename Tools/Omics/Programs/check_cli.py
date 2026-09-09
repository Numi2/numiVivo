#!/usr/bin/env python3
"""Program-scoring publication, replay, resident/streamed equivalence and rejections."""
import argparse,copy,hashlib,json,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--prepared',type=Path,required=True)
p.add_argument('--programs',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--imported',type=Path)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def write(p,x): p.write_text(json.dumps(x,indent=2)+'\n')
commands=[]
def run(*args,fail=None):
    start=time.monotonic();r=subprocess.run([str(a.binary.resolve()),*map(str,args)],capture_output=True,text=True,timeout=600)
    i=len(commands);(a.out/f'{i:02d}.stdout').write_text(r.stdout);(a.out/f'{i:02d}.stderr').write_text(r.stderr)
    commands.append(dict(arguments=list(map(str,args)),exitCode=r.returncode,seconds=time.monotonic()-start,expectedFailure=fail))
    write(a.out/'commands.json',commands)
    if fail is None: assert r.returncode==0,r.stderr[-4000:]
    else: assert r.returncode==65 and fail in r.stderr.lower(),r.stderr[-4000:]
    return r
mapping=json.loads((a.prepared/'mapping.json').read_text())
programs=json.loads(a.programs.read_text())
plan=dict(schemaVersion=1,mapping=mapping,contrasts=[],programs=programs,reduction={})
write(a.out/'plan.json',plan)
source=a.prepared/'prepared.h5ad'
run('singlecell-h5ad-pseudobulk',source,'--plan',a.out/'plan.json','--output',a.out/'bundle')
run('singlecell-h5ad-pseudobulk-verify',a.out/'bundle')
run('singlecell-h5ad-pseudobulk',source,'--plan',a.out/'plan.json','--output',a.out/'repeat')
assert (a.out/'bundle/receipt.json').read_bytes()==(a.out/'repeat/receipt.json').read_bytes()
checks=['streamed programs with PCA','native reconstruction','exact repeat']
if a.imported:
    store=a.out/'store';counts=a.out/'counts.json';receipt=a.out/'analysis.json'
    run('singlecell-run',a.imported/'manifest.json','--store',store,'--output',counts)
    resident=dict(schemaVersion=1,id='program-scoring-qualification',normalizationTarget=10000,contrasts=[],programs=programs)
    write(a.out/'resident-plan.json',resident)
    run('singlecell-analyze',counts,'--plan',a.out/'resident-plan.json','--store',store,'--output',receipt)
    run('singlecell-analysis-verify',receipt,'--store',store)
    run('singlecell-analysis-export',receipt,'--store',store,'--output',a.out/'resident-report.json')
    r=json.loads((a.out/'resident-report.json').read_text())['programs'];s=json.loads((a.out/'bundle/report.json').read_text())['programs']
    def rows(x):return {(c['sampleID'],c['barcode']):(v,n) for c,v,n in zip(x['cells'],x['scores'],x['detectedMembers'],strict=True)}
    assert rows(r)==rows(s)
    for k in ['programs','method','normalizationTarget','updates']:assert r[k]==s[k],k
    checks.append('resident and streamed scores/detection exactly equal by cell identity')
# All are required controlled failures, with no output bundle/staging left behind.
for case in ['missing','organism','score-budget','update-budget','duplicate','provenance']:
    bad=copy.deepcopy(plan);d=bad['programs']['definitions'][0]
    if case=='missing':d['members'].append(dict(featureID='NUMIVIVO_UNMATCHED_MARKER',weight=1));d['minimumWeightCoverage']=1
    elif case=='organism':d['organism']='NCBITaxon:10090'
    elif case=='score-budget':bad['programs']['maximumScoreValues']=1
    elif case=='update-budget':bad['programs']['maximumUpdates']=1
    elif case=='duplicate':d['members'].append(copy.deepcopy(d['members'][0]))
    elif case=='provenance':d['sourceVersion']=''
    path=a.out/f'bad-{case}.json';write(path,bad)
    expected={'missing':'missing-feature','organism':'organism','score-budget':'score-value budget','update-budget':'sparse-update budget','duplicate':'program definition','provenance':'program definition'}[case]
    run('singlecell-h5ad-pseudobulk',source,'--plan',path,'--output',a.out/f'rejected-{case}',fail=expected)
    assert not (a.out/f'rejected-{case}').exists()
    checks.append('controlled '+case+' rejection')
assert not list(a.out.glob('.numivivo-stream-*'))
checks.append('failed publication leaves no staging bundle')
write(a.out/'checks.json',dict(status='passed',checks=checks,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),qualification='Software scoring and replay; biological observations are evaluated separately'))
print(json.dumps(dict(status='passed',checks=checks),indent=2))
