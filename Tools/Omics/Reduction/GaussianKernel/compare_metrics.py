#!/usr/bin/env python3
"""Compare every original biological metric and complete neighbor/prediction array."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--original',type=Path,required=True);a=p.parse_args()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
results=[]
for c in ('hagai','kang','ding'):
 r=a.root/c;old=a.original/('evaluation-'+c);d=json.loads((r/'checks.json').read_text());changes=[]
 def walk(x,y,path=''):
  if isinstance(x,dict):
   assert set(x)==set(y)
   for k in x:
    if k not in ('evaluationSeconds','seconds'):walk(x[k],y[k],path+'/'+k)
  elif isinstance(x,list):
   assert len(x)==len(y)
   for i,(v,w) in enumerate(zip(x,y)):walk(v,w,path+'/'+str(i))
  elif x!=y:changes.append(dict(path=path,original=y,candidate=x))
 walk(d['metrics'],d['originalNativeMetrics']);assert d['gates']==d['originalNativeGates'];arrays=[]
 for n in (['neighbors','predictions'] if c=='ding' else ['neighbors']):
  xp=r/('candidate-'+n+'.npz');yp=old/('native-'+n+'.npz');x=np.load(xp,allow_pickle=False);y=np.load(yp,allow_pickle=False);assert set(x.files)==set(y.files)
  for k in x.files:
   assert x[k].shape==y[k].shape;changed=int(np.count_nonzero(x[k]!=y[k]));arrays.append(dict(file=n,key=k,entries=x[k].size,changedEntries=changed))
 results.append(dict(cohort=c,allGatesIdentical=True,metricChanges=changes,arrays=arrays,candidateChecksSHA256=sha(r/'checks.json'),originalChecksSHA256=sha(old/'checks.json')))
out=dict(status='compared-all-values',cohorts=results,checkerSHA256=sha(Path(__file__)))
(a.root/'original-comparison.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
