#!/usr/bin/env python3
"""Full real integrated PCA -> HNSW -> clustering / embedding lifecycle."""
import argparse,json,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['binary','input','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def save(path,v):path.write_text(json.dumps(v,indent=2)+'\n')
def run(label,args):
 command=['/usr/bin/time','-l',str(a.binary),*map(str,args)];r=subprocess.run(command,capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
 rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);seconds=re.search(r'([\d.]+) real',r.stderr)
 commands.append(dict(label=label,command=command,exitCode=r.returncode,seconds=float(seconds[1]) if seconds else None,maximumResidentBytes=int(rss[1]) if rss else None));save(a.out/'commands.json',commands);assert r.returncode==0,r.stderr;print(label+' passed',flush=True)
plan=dict(schemaVersion=1,inputKind='integrated',neighbors=dict(neighbors=15,representation='integrated'),execution=dict(workers=1),approximation=dict(connections=16,constructionWidth=200,searchWidth=128,seed=7,maximumDistanceEvaluations=500000000,scoreCacheBytes=33554432),storage='binary')
save(a.out/'graph-plan.json',plan);run('graph',['singlecell-pca-neighbors',a.input,'--plan',a.out/'graph-plan.json','--output',a.out/'graph']);run('graph-verify',['singlecell-pca-neighbors-verify',a.out/'graph'])
for action,options in [('cluster',dict(clustering=dict())),('embed',dict(embedding=dict(epochs=500,seed=7,maximumUpdates=1000000000)))]:
 path=a.out/(action+'-plan.json');save(path,dict(schemaVersion=1,**options));run(action,['singlecell-graph-'+action,a.out/'graph','--plan',path,'--output',a.out/action]);run(action+'-verify',['singlecell-graph-'+action+'-verify',a.out/action])
save(a.out/'checks.json',dict(status='passed',commands=len(commands),representation='integrated',qualification='Full-cohort integrated graph, clustering and embedding lifecycle; HNSW recall and biological preservation require separate evaluation.'))
