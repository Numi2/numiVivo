#!/usr/bin/env python3
"""Run complete real-data fitted/query graph lifecycles; stdlib, native CPU work only."""
import argparse,hashlib,json,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[];results=[]
def run(label,args):
    command=['/usr/bin/time','-l',str(a.binary),*map(str,args)]
    r=subprocess.run(command,capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
    rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);elapsed=re.search(r'([\d.]+) real',r.stderr)
    commands.append(dict(label=label,command=command,exitCode=r.returncode,seconds=float(elapsed[1]) if elapsed else None,maximumResidentBytes=int(rss[1]) if rss else None))
    (a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');assert r.returncode==0,r.stderr
for case in json.loads(a.inputs.read_text()):
    name=case['name'];root=a.out/name;root.mkdir();source=Path(case['input'])
    if case['kind']=='fitted':
        run(name+'-fit',['singlecell-h5ad-pca',source,'--plan',case['fitPlan'],'--output',root/'pca']);source=root/'pca'
    elif 'trainingInput' in case:
        run(name+'-query-fit',['singlecell-h5ad-pca',case['trainingInput'],'--plan',case['fitPlan'],'--output',root/'reference'])
        run(name+'-query-project',['singlecell-h5ad-pca-query',source,'--plan',case['queryPlan'],'--reference',root/'reference','--output',root/'query'])
        source=root/'query'
    for workers in [1,4]:
        plan=dict(schemaVersion=1,inputKind=case['kind'],neighbors=dict(neighbors=15,maximumDistancePairs=200_000_000),execution=dict(workers=workers,queryBlockRows=128,candidateBlockRows=2048))
        path=root/('plan-'+str(workers)+'.json');path.write_text(json.dumps(plan,indent=2)+'\n')
        run(name+'-workers-'+str(workers),['singlecell-pca-neighbors',source,'--plan',path,'--output',root/('workers-'+str(workers))])
    serial=(root/'workers-1/graph.json').read_bytes();parallel=(root/'workers-4/graph.json').read_bytes();assert serial==parallel
    run(name+'-verify',['singlecell-pca-neighbors-verify',root/'workers-4'])
    graph=json.loads(parallel);results.append(dict(name=name,inputKind=case['kind'],cells=len(graph['cells']),graphSHA256=hashlib.sha256(parallel).hexdigest(),serialParallelGraphBytesExact=True,uniqueDistancePairs=graph['distancePairs'],directedDistanceEvaluations=graph['distancePairs']*2))
    print(json.dumps(results[-1]),flush=True)
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',cases=results,commands=len(commands)),indent=2)+'\n')
