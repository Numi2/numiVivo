#!/usr/bin/env python3
"""Native full-cohort HNSW publication/reconstruction; reference recall is separate."""
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
    print(label+' passed',flush=True)
for case in json.loads(a.inputs.read_text()):
    name=case['name'];root=a.out/name;root.mkdir()
    run(name+'-fit',['singlecell-h5ad-pca',case['source'],'--plan',case['fitPlan'],'--output',root/'pca'])
    plan=dict(schemaVersion=1,inputKind='fitted',neighbors=dict(neighbors=15),execution=dict(workers=1),approximation=dict(connections=16,constructionWidth=200,searchWidth=128,seed=7,maximumDistanceEvaluations=500_000_000,scoreCacheBytes=33_554_432))
    path=root/'plan.json';path.write_text(json.dumps(plan,indent=2)+'\n')
    run(name+'-hnsw',['singlecell-pca-neighbors',root/'pca','--plan',path,'--output',root/'graph'])
    run(name+'-verify',['singlecell-pca-neighbors-verify',root/'graph'])
    graph=(root/'graph/graph.json').read_bytes();execution=json.loads((root/'graph/execution.json').read_text())
    result=dict(name=name,cells=execution['cells'],graphSHA256=hashlib.sha256(graph).hexdigest(),distanceEvaluations=execution['directedDistanceEvaluations'],hnsw=execution['hnsw'])
    results.append(result);print(json.dumps(result),flush=True)
    (a.out/'completed.json').write_text(json.dumps(results,indent=2)+'\n')
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',cases=results,commands=len(commands),scope='Native lifecycle only; independent recall checks required'),indent=2)+'\n')
