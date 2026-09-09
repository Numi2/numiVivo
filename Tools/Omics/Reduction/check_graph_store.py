#!/usr/bin/env python3
"""Fitted/query and exact/approximate graph stores, reconstruction and tampering."""
import argparse,hashlib,json,shutil,struct,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--fitted',type=Path,required=True);p.add_argument('--query',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def run(label,args,ok=True):
    command=[str(a.binary),*map(str,args)];r=subprocess.run(command,capture_output=True,text=True)
    (a.out/(label+'.log')).write_text(r.stdout+r.stderr);commands.append(dict(label=label,command=command,exitCode=r.returncode,expectedSuccess=ok))
    (a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');assert (r.returncode==0)==ok,(label,r.stderr)
def compare(root,base):
    report=json.loads((root/'graph.json').read_text());g=json.loads((base/'graph.json').read_text());n=len(g['cells']);k=g['options']['neighbors']
    def records(name):return list(struct.iter_unpack('<IId',(root/name).read_bytes()))
    assert records('neighbors.bin')==[(i//k,g['neighborIndices'][i],v) for i,v in enumerate(g['neighborDistances'])]
    assert records('bandwidths.bin')==[(i,j,g[key][i]) for i in range(n) for j,key in enumerate(['rhos','sigmas','kernelMassResiduals'])]
    assert list(struct.iter_unpack('<IIQ',(root/'offsets.bin').read_bytes()))==[(i,0,v) for i,v in enumerate(g['rowOffsets'])]
    assert records('edges.bin')==[(i,g['columnIndices'][j],g['weights'][j]) for i in range(n) for j in range(g['rowOffsets'][i],g['rowOffsets'][i+1])]
    for key in ['connectedComponents','isolatedCells','distancePairs','options','dimensions','method']:assert report[key]==g[key],key
for kind,input in [('fitted',a.fitted),('query',a.query)]:
    for mode in ['exact','hnsw']:
        name=kind+'-'+mode;root=a.out/name;root.mkdir()
        plan=dict(schemaVersion=1,inputKind=kind,neighbors=dict(neighbors=5),execution=dict(workers=2 if mode=='exact' else 1),storage='binary')
        if mode=='hnsw':plan['approximation']=dict(connections=8,constructionWidth=64,searchWidth=64)
        path=root/'plan.json';path.write_bytes(canonical(plan))
        run(name+'-publish',['singlecell-pca-neighbors',input,'--plan',path,'--output',root/'binary'])
        run(name+'-verify',['singlecell-pca-neighbors-verify',root/'binary'])
        run(name+'-repeat',['singlecell-pca-neighbors',input,'--plan',path,'--output',root/'repeat'])
        assert (root/'binary/receipt.json').read_bytes()==(root/'repeat/receipt.json').read_bytes()
        del plan['storage'];jsonplan=root/'json-plan.json';jsonplan.write_bytes(canonical(plan))
        run(name+'-json',['singlecell-pca-neighbors',input,'--plan',jsonplan,'--output',root/'json'])
        compare(root/'binary',root/'json')
base=a.out/'fitted-hnsw/binary'
for key in ['neighbors','bandwidths','offsets','edges']:
    root=a.out/('tamper-'+key);shutil.copytree(base,root)
    path=root/(key+'.bin');data=bytearray(path.read_bytes());data[8]^=1;path.write_bytes(data)
    graph=json.loads((root/'graph.json').read_text());graph[key]=dict(bytes=list(hashlib.sha256(data).digest()));(root/'graph.json').write_bytes(canonical(graph))
    receipt=json.loads((root/'receipt.json').read_text());receipt['graph']=dict(bytes=list(hashlib.sha256((root/'graph.json').read_bytes()).digest()));(root/'receipt.json').write_bytes(canonical(receipt))
    run('rehashed-'+key,['singlecell-pca-neighbors-verify',root],False)
run('overwrite',['singlecell-pca-neighbors',a.fitted,'--plan',a.out/'fitted-hnsw/plan.json','--output',base],False)
bad=json.loads((a.out/'fitted-hnsw/plan.json').read_text());bad['storage']='unknown';path=a.out/'invalid-storage.json';path.write_bytes(canonical(bad))
run('invalid-storage',['singlecell-pca-neighbors',a.fitted,'--plan',path,'--output',a.out/'invalid-output'],False)
assert not any(p.name.startswith('.numivivo-pca-neighbors-') for p in a.out.rglob('*'))
result=dict(status='passed',commands=len(commands),expectedRejections=sum(not c['expectedSuccess'] for c in commands),fittedAndQuery=True,exactAndHNSW=True,allBinaryNumericsExact=True,replayExact=True,allFourRehashedFilesRejected=True,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest())
(a.out/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
