#!/usr/bin/env python3
"""Full-cohort file clustering lifecycle and frozen pre-refactor numerical oracle."""
import argparse,json,re,subprocess,hashlib
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--oracle',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[];results=[]
clustering=dict(resolution=1,seed=7,maximumSweeps=100,maximumLevels=32,levelTolerance=1e-7)
(a.out/'clustering-plan.json').write_text(json.dumps(dict(schemaVersion=1,clustering=clustering,maximumEdgeVisits=1000000000),indent=2)+'\n')
(a.out/'legacy-options.json').write_text(json.dumps(clustering,indent=2)+'\n')
def run(label,binary,args):
    command=['/usr/bin/time','-l',str(binary),*map(str,args)];r=subprocess.run(command,capture_output=True,text=True)
    (a.out/(label+'.log')).write_text(r.stdout+r.stderr);rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);seconds=re.search(r'([\d.]+) real',r.stderr)
    commands.append(dict(label=label,command=command,exitCode=r.returncode,seconds=float(seconds[1]) if seconds else None,maximumResidentBytes=int(rss[1]) if rss else None))
    (a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');assert r.returncode==0,r.stderr;print(label+' passed',flush=True)
for case in json.loads(a.inputs.read_text()):
    name=case['name'];root=a.out/name;root.mkdir()
    run(name+'-fit',a.binary,['singlecell-h5ad-pca',case['source'],'--plan',case['fitPlan'],'--output',root/'pca'])
    plan=dict(schemaVersion=1,inputKind='fitted',neighbors=dict(neighbors=15),execution=dict(workers=1),approximation=dict(connections=16,constructionWidth=200,searchWidth=128,seed=7,maximumDistanceEvaluations=500000000,scoreCacheBytes=33554432),storage='binary')
    path=root/'graph-plan.json';path.write_text(json.dumps(plan,indent=2)+'\n')
    run(name+'-graph',a.binary,['singlecell-pca-neighbors',root/'pca','--plan',path,'--output',root/'graph'])
    run(name+'-cluster',a.binary,['singlecell-graph-cluster',root/'graph','--plan',a.out/'clustering-plan.json','--output',root/'cluster'])
    run(name+'-verify',a.binary,['singlecell-graph-cluster-verify',root/'cluster'])
    run(name+'-legacy-oracle',a.oracle,[case['legacyGraph'],a.out/'legacy-options.json',root/'legacy-result.json'])
    raw=(root/'cluster/result.json').read_bytes();assert raw==(root/'legacy-result.json').read_bytes(),'Frozen legacy result differs; retain evidence before investigating'
    value=json.loads(raw);entry=dict(case=name,cells=len(value['cells']),clusters=len(value['clusterSizes']),modularity=value['modularity'],levels=value['levels'],disconnectedCommunities=value['disconnectedCommunities'],allLegacyResultBytesExact=True,resultSHA256=hashlib.sha256(raw).hexdigest())
    results.append(entry);(a.out/'completed.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(entry),flush=True)
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',commands=len(commands),cases=results,scope='Native lifecycle and exact frozen algorithm comparison; independent graph/partition qualification remains separate'),indent=2)+'\n')
