#!/usr/bin/env python3
"""Full-cohort binary graph lifecycle with a same-executable JSON baseline."""
import argparse,json,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def run(label,args):
    command=['/usr/bin/time','-l',str(a.binary),*map(str,args)]
    r=subprocess.run(command,capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
    rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);seconds=re.search(r'([\d.]+) real',r.stderr)
    commands.append(dict(label=label,command=command,exitCode=r.returncode,seconds=float(seconds[1]) if seconds else None,maximumResidentBytes=int(rss[1]) if rss else None))
    (a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');assert r.returncode==0,r.stderr
    print(label+' passed',flush=True)
for case in json.loads(a.inputs.read_text()):
    name=case['name'];root=a.out/name;root.mkdir()
    run(name+'-fit',['singlecell-h5ad-pca',case['source'],'--plan',case['fitPlan'],'--output',root/'pca'])
    plan=dict(schemaVersion=1,inputKind='fitted',neighbors=dict(neighbors=15),execution=dict(workers=1),approximation=dict(connections=16,constructionWidth=200,searchWidth=128,seed=7,maximumDistanceEvaluations=500000000,scoreCacheBytes=33554432),storage='binary')
    path=root/'binary-plan.json';path.write_text(json.dumps(plan,indent=2)+'\n')
    run(name+'-binary',['singlecell-pca-neighbors',root/'pca','--plan',path,'--output',root/'binary'])
    run(name+'-verify',['singlecell-pca-neighbors-verify',root/'binary'])
    del plan['storage'];path=root/'json-plan.json';path.write_text(json.dumps(plan,indent=2)+'\n')
    run(name+'-json',['singlecell-pca-neighbors',root/'pca','--plan',path,'--output',root/'json'])
    if name=='baron':
        del plan['approximation'];plan['storage']='binary';plan['execution']['workers']=4;plan['neighbors']['maximumDistancePairs']=200000000
        path=root/'exact-plan.json';path.write_text(json.dumps(plan,indent=2)+'\n')
        run(name+'-exact-binary',['singlecell-pca-neighbors',root/'pca','--plan',path,'--output',root/'exact-binary'])
        run(name+'-exact-verify',['singlecell-pca-neighbors-verify',root/'exact-binary'])
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',commands=len(commands),scope='Native lifecycle only; independent binary graph and baseline checks required'),indent=2)+'\n')
