#!/usr/bin/env python3
"""Check that one global PCA scale preserves frozen baseline readouts."""
import argparse,json
from pathlib import Path
import numpy as np
from check_full_integration import metrics as donor_metrics
from evaluate_ding_integration import metrics as ding_metrics
p=argparse.ArgumentParser();p.add_argument('--cohort',required=True,choices=['kang','hagai','ding']);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
if a.cohort=='ding':
 inp=np.load(a.inputs/'inputs.npz',allow_pickle=False);x=inp['scores'];baseline=json.loads((a.inputs/'baseline.json').read_text());old_neighbors=a.inputs/'baseline-neighbors.npz'
else:
 inp=np.load(a.inputs/'evaluation-inputs.npz',allow_pickle=False);x=inp['scores'];baseline=json.loads((a.inputs/'checks.json').read_text())['baseline'];old_neighbors=a.inputs/'baseline-neighbors.npz'
scale=float(np.median(np.linalg.norm(x,axis=1)));assert scale>0
if a.cohort=='ding':m=ding_metrics(x/scale,inp,30,4,a.out,'baseline')
else:m=donor_metrics(x/scale,inp['donors'],inp['conditions'],inp['cellTypes'],inp['program'],30,a.out/'baseline-neighbors.npz')
(a.out/'baseline.json').write_text(json.dumps(m,indent=2,allow_nan=False)+'\n')
errors=[]
def compare(a,b,path=''):
 if isinstance(a,dict):
  for key in a.keys()&b.keys():
   if key not in ['evaluationSeconds']:compare(a[key],b[key],path+'/'+key)
 elif isinstance(a,list):
  assert len(a)==len(b),path
  for i,(v,w) in enumerate(zip(a,b)):compare(v,w,path+'/'+str(i))
 elif isinstance(a,(float,int)):
  if abs(a-b)>1e-12:errors.append(dict(path=path,baseline=a,scaled=b,absoluteError=abs(a-b)))
 else:assert a==b,path
compare(baseline,m)
old=np.load(old_neighbors);new=np.load(a.out/'baseline-neighbors.npz');neighbor_changes={key:int(np.sum(old[key]!=new[key])) for key in old.files};assert all(value==0 for value in neighbor_changes.values())
assert not errors,errors
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',cohort=a.cohort,cells=len(x),globalScale=scale,allBaselineMetricDifferencesAtMost1eMinus12=True,neighborIndexChanges=neighbor_changes,qualification='Uniform positive scaling leaves every retained exact neighbor index and original baseline metric unchanged within 1e-12. Timing excluded; all source rows retained.'),indent=2)+'\n')
