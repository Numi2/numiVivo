#!/usr/bin/env python3
"""HNSW fitted/query lifecycle, repeatability, cache invariance and negative controls."""
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
    plan=dict(schemaVersion=1,inputKind=kind,neighbors=dict(neighbors=5),execution=dict(workers=1),approximation=dict(connections=8,constructionWidth=64,searchWidth=64))
    path=a.out/(kind+'.json');write(path,plan);out=a.out/kind
    run(['singlecell-pca-neighbors',source,'--plan',path,'--output',out]);run(['singlecell-pca-neighbors-verify',out])
    run(['singlecell-pca-neighbors',source,'--plan',path,'--output',a.out/(kind+'-repeat')]);assert (out/'receipt.json').read_bytes()==(a.out/(kind+'-repeat/receipt.json')).read_bytes()
    cached=copy.deepcopy(plan);cached['approximation']['scoreCacheBytes']=262144;write(a.out/(kind+'-cache.json'),cached)
    run(['singlecell-pca-neighbors',source,'--plan',a.out/(kind+'-cache.json'),'--output',a.out/(kind+'-cache')]);assert (out/'graph.json').read_bytes()==(a.out/(kind+'-cache/graph.json')).read_bytes()
    # Large search width covers this entire fixture; compare the unchanged exact route.
    exact=copy.deepcopy(plan);exact.pop('approximation');write(a.out/(kind+'-exact.json'),exact)
    run(['singlecell-pca-neighbors',source,'--plan',a.out/(kind+'-exact.json'),'--output',a.out/(kind+'-exact')])
    graph=json.loads((out/'graph.json').read_text());baseline=json.loads((a.out/(kind+'-exact/graph.json')).read_text())
    for key in ['neighborIndices','neighborDistances','rhos','sigmas','kernelMassResiduals','rowOffsets','columnIndices','weights','connectedComponents','isolatedCells']:assert graph[key]==baseline[key],(kind,key)
    run(['singlecell-pca-neighbors',source,'--plan',path,'--output',out],'destination exists')
plan=json.loads((a.out/'fitted.json').read_text())
for name,edit,error in [('workers',lambda v:v['execution'].update(workers=2),'one serial worker'),('cache',lambda v:v['approximation'].update(scoreCacheBytes=1),'resource bounds'),('search',lambda v:v['approximation'].update(searchWidth=2),'resource bounds'),('degree',lambda v:v['approximation'].update(connections=65),'resource bounds'),('budget',lambda v:v['approximation'].update(maximumDistanceEvaluations=1),'distance-evaluation budget')]:
    bad=copy.deepcopy(plan);edit(bad);path=a.out/(name+'.json');write(path,bad)
    run(['singlecell-pca-neighbors',a.fitted,'--plan',path,'--output',a.out/('reject-'+name)],error);assert not (a.out/('reject-'+name)).exists()
for kind in ['graph','execution','input']:
    dst=a.out/('tamper-'+kind);shutil.copytree(a.out/'fitted',dst)
    if kind=='input':
        path=dst/'input/scores.bin';value=bytearray(path.read_bytes());value[0]^=1;path.write_bytes(value);receiptPath=dst/'input/receipt.json';key='scores';error='PCA source reconstruction differs'
    else:
        path=dst/(kind+'.json');value=json.loads(path.read_text())
        if kind=='graph':value['weights'][0]*=0.5
        else:value['hnsw']['constructionDistances']+=1
        write(path,value);receiptPath=dst/'receipt.json';key='graph' if kind=='graph' else 'executionReport';error='PCA neighbor source reconstruction differs'
    receipt=json.loads(receiptPath.read_text());receipt[key]={'bytes':list(hashlib.sha256(path.read_bytes()).digest())};write(receiptPath,receipt)
    run(['singlecell-pca-neighbors-verify',dst],error)
assert not list(a.out.glob('.numivivo-pca-neighbors-*'))
write(a.out/'checks.json',dict(status='passed',commands=len(commands),expectedRejections=sum(c['expectedRejection'] is not None for c in commands),fittedAndQueryInputs=True,replayExact=True,cacheInvariantGraphBytes=True,smallFixtureExactGraphNumerics=True,rehashedTamperingRejected=True,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest()))
