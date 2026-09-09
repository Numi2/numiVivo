#!/usr/bin/env python3
"""Fitted/query, exact/HNSW embedding lifecycle, frozen oracle and rejection controls."""
import argparse,hashlib,json,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--oracle',type=Path,required=True);p.add_argument('--graph-fixtures',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':')).encode()
def run(label,binary,args,ok=True):
    command=[str(binary),*map(str,args)];r=subprocess.run(command,capture_output=True,text=True)
    (a.out/(label+'.log')).write_text(r.stdout+r.stderr);commands.append(dict(label=label,command=command,exitCode=r.returncode,expectedSuccess=ok))
    (a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');assert (r.returncode==0)==ok,(label,r.stderr)
plan=dict(schemaVersion=1,embedding=dict(dimensions=2,epochs=20,minimumDistance=0.1,spread=1,learningRate=1,negativeSampleRate=5,repulsionStrength=1,seed=7,maximumUpdates=1000000))
path=a.out/'plan.json';path.write_bytes(canonical(plan));opts=a.out/'legacy-options.json';opts.write_bytes(canonical(plan['embedding']))
for kind in ['fitted','query']:
    for mode in ['exact','hnsw']:
        name=kind+'-'+mode;source=a.graph_fixtures/name
        plan['embedding']['dimensions']=2 if mode=='exact' else 3
        path.write_bytes(canonical(plan));opts.write_bytes(canonical(plan['embedding']))
        run(name+'-publish',a.binary,['singlecell-graph-embed',source/'binary','--plan',path,'--output',a.out/name])
        run(name+'-verify',a.binary,['singlecell-graph-embed-verify',a.out/name])
        run(name+'-repeat',a.binary,['singlecell-graph-embed',source/'binary','--plan',path,'--output',a.out/(name+'-repeat')])
        assert (a.out/name/'receipt.json').read_bytes()==(a.out/(name+'-repeat')/'receipt.json').read_bytes()
        run(name+'-legacy',a.oracle,[source/'json/graph.json',source/'json/input/scores.bin',opts,a.out/(name+'-legacy.json')])
        assert (a.out/name/'result.json').read_bytes()==(a.out/(name+'-legacy.json')).read_bytes()
base=a.out/'fitted-hnsw'
for name,key,field in [('result.json','result','coordinates'),('execution.json','executionReport','scheduleBytes')]:
    root=a.out/('tamper-'+key);shutil.copytree(base,root);value=json.loads((root/name).read_text())
    if field=='coordinates':value[field][0][0]+=1
    else:value[field]+=1
    (root/name).write_bytes(canonical(value));receipt=json.loads((root/'receipt.json').read_text());receipt[key]=dict(bytes=list(hashlib.sha256((root/name).read_bytes()).digest()));(root/'receipt.json').write_bytes(canonical(receipt))
    run('rehashed-'+key,a.binary,['singlecell-graph-embed-verify',root],False)
root=a.out/'tamper-parent';shutil.copytree(base,root);edge=root/'input/edges.bin';data=bytearray(edge.read_bytes());data[8]^=1;edge.write_bytes(data)
graph=json.loads((root/'input/graph.json').read_text());graph['edges']=dict(bytes=list(hashlib.sha256(data).digest()));(root/'input/graph.json').write_bytes(canonical(graph))
gr=json.loads((root/'input/receipt.json').read_text());gr['graph']=dict(bytes=list(hashlib.sha256((root/'input/graph.json').read_bytes()).digest()));(root/'input/receipt.json').write_bytes(canonical(gr))
receipt=json.loads((root/'receipt.json').read_text());receipt['input']=dict(bytes=list(hashlib.sha256((root/'input/receipt.json').read_bytes()).digest()));(root/'receipt.json').write_bytes(canonical(receipt))
run('rehashed-parent',a.binary,['singlecell-graph-embed-verify',root],False)
low=json.loads(json.dumps(plan));low['embedding']['maximumUpdates']=1;lowpath=a.out/'low-work.json';lowpath.write_bytes(canonical(low))
run('low-work',a.binary,['singlecell-graph-embed',a.graph_fixtures/'fitted-hnsw/binary','--plan',lowpath,'--output',a.out/'rejected-work'],False)
run('resident-input',a.binary,['singlecell-graph-embed',a.graph_fixtures/'fitted-hnsw/json','--plan',path,'--output',a.out/'rejected-storage'],False)
run('overwrite',a.binary,['singlecell-graph-embed',a.graph_fixtures/'fitted-hnsw/binary','--plan',path,'--output',base],False)
assert not any(p.name.startswith(('.embedding-schedule-','.numivivo-graph-embedding-')) for p in a.out.rglob('*'))
result=dict(status='passed',commands=len(commands),expectedRejections=sum(not v['expectedSuccess'] for v in commands),fittedAndQuery=True,exactAndHNSW=True,twoAndThreeDimensions=True,allFrozenOracleBytesExact=True,replayExact=True,rehashedResultWorkAndParentRejected=True,failedStagingRemoved=True,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),oracleSHA256=hashlib.sha256(a.oracle.read_bytes()).hexdigest())
(a.out/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
