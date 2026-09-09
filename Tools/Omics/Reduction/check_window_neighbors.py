#!/usr/bin/env python3
"""File-backed fitted/query graph lifecycle, tile invariance and tamper controls."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--fitted',type=Path,required=True);p.add_argument('--query',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def write(path,value):path.write_text(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False))
def run(args,error=None):
    r=subprocess.run([str(a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(str(len(commands))+'.log')).write_text(r.stdout+r.stderr)
    commands.append(dict(arguments=list(map(str,args)),exitCode=r.returncode,expectedRejection=error));write(a.out/'commands.json',commands)
    assert (r.returncode!=0 and error in r.stderr) if error else r.returncode==0,r.stderr
for kind,source in [('fitted',a.fitted),('query',a.query)]:
    baseline=None
    for workers,qrows,crows in [(1,3,5),(2,7,13),(4,9,32)]:
        plan=dict(schemaVersion=1,inputKind=kind,neighbors=dict(neighbors=5),execution=dict(workers=workers,queryBlockRows=qrows,candidateBlockRows=crows))
        path=a.out/(kind+'-'+str(workers)+'.json');write(path,plan);out=a.out/(kind+'-'+str(workers))
        run(['singlecell-pca-neighbors',source,'--plan',path,'--output',out])
        graph=(out/'graph.json').read_bytes()
        if baseline is None:baseline=graph
        else:assert baseline==graph
    run(['singlecell-pca-neighbors-verify',out])
    run(['singlecell-pca-neighbors',source,'--plan',path,'--output',a.out/(kind+'-repeat')])
    assert (out/'receipt.json').read_bytes()==(a.out/(kind+'-repeat/receipt.json')).read_bytes()
    run(['singlecell-pca-neighbors',source,'--plan',path,'--output',out],'destination exists')
plan=dict(schemaVersion=1,inputKind='fitted',neighbors=dict(neighbors=5),execution=dict(workers=4))
for name,edit,error in [('workers',lambda v:v['execution'].update(workers=0),'worker or tile bounds'),('tile',lambda v:v['execution'].update(candidateBlockRows=8193),'worker or tile bounds'),('budget',lambda v:v['neighbors'].update(maximumDistancePairs=1),'distance-pair budget'),('representation',lambda v:v['neighbors'].update(representation='integrated'),'schema or representation')]:
    bad=copy.deepcopy(plan);edit(bad);path=a.out/(name+'.json');write(path,bad)
    run(['singlecell-pca-neighbors',a.fitted,'--plan',path,'--output',a.out/('reject-'+name)],error);assert not (a.out/('reject-'+name)).exists()
for kind in ['graph','execution','input']:
    dst=a.out/('tamper-'+kind);shutil.copytree(a.out/'fitted-4',dst)
    if kind=='input':
        path=dst/'input/scores.bin';value=bytearray(path.read_bytes());value[8]^=1;path.write_bytes(value);receiptPath=dst/'input/receipt.json';key='scores';error='PCA source reconstruction differs'
    else:
        path=dst/(kind+'.json');value=json.loads(path.read_text())
        if kind=='graph':value['weights'][0]*=0.5
        else:value['directedDistanceEvaluations']+=1
        write(path,value);receiptPath=dst/'receipt.json';key='graph' if kind=='graph' else 'executionReport';error='PCA neighbor source reconstruction differs'
    receipt=json.loads(receiptPath.read_text());receipt[key]={'bytes':list(hashlib.sha256(path.read_bytes()).digest())};write(receiptPath,receipt)
    run(['singlecell-pca-neighbors-verify',dst],error)
assert not list(a.out.glob('.numivivo-pca-neighbors-*'))
write(a.out/'checks.json',dict(status='passed',commands=len(commands),expectedRejections=sum(c['expectedRejection'] is not None for c in commands),fittedAndQueryInputs=True,allSerialParallelGraphsBytesExact=True,replayExact=True,rehashedTamperingRejected=True,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest()))
