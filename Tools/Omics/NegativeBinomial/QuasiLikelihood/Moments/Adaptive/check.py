#!/usr/bin/env python3
"""Independently check adaptive outputs with bounded direct PMF summation."""
import argparse,gzip,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from scipy.special import gammaln,xlogy
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--binary',required=True);p.add_argument('--old',type=Path,required=True);a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True)
sha=lambda b:hashlib.sha256(b).hexdigest()
def check(mu,alpha,m):
 mass=first=second=0.;r=1/alpha
 for start in range(m['firstCount'],m['lastCount']+1,250000):
  k=np.arange(start,min(start+250000,m['lastCount']+1),dtype=float)
  prob=np.exp(gammaln(k+r)-gammaln(k+1)-gammaln(r)-r*np.log1p(mu/r)-k*np.log1p(r/mu))
  d=2*(xlogy(k,k/mu)-(k+r)*np.log1p((k-mu)/(mu+r)))
  if start==0:prob[0]=np.exp(-r*np.log1p(mu/r));d[0]=2*r*np.log1p(mu/r)
  mass+=float(prob.sum());first+=float(prob@d);second+=float(prob@(d*d))
 mean=first/mass;variance=second/mass-mean*mean
 expected=dict(mean=mean,variance=variance,devianceScale=2*mean/variance,degreesOfFreedom=2*mean*mean/variance)
 errors={k:abs(m[k]/v-1) for k,v in expected.items()}
 b=m.get('summation',{});bound=max((m['meanTruncationBound']+b.get('meanErrorBound',0))/m['mean'],(m['varianceTruncationBound']+b.get('varianceErrorBound',0))/m['variance'])
 return dict(reference=expected,relativeErrors=errors,referenceMass=mass,combinedRelativeBound=bound,evaluatedCounts=m['evaluatedCounts'],coveredCounts=m['lastCount']-m['firstCount']+1,status='passed' if max(errors.values())<=2e-7 and abs(mass-1)<=1e-7+m['omittedProbabilityBound'] and bound<=1e-10 else 'failed')
oldInputs=json.loads((a.old/'work-sensitivity/input.json').read_text());oldItems=json.loads((a.old/'work-sensitivity/protocol.json').read_text())['items'];oldOutputs=json.loads(gzip.decompress((a.old/'work-sensitivity/native.json.gz').read_bytes()))
native=json.loads((a.root/'twenty-native.json').read_text())['residuals'];assert len(native)==len(oldInputs)==20
checks=[];failures=[]
for item,i,v,old in zip(oldItems,oldInputs,native,oldOutputs,strict=True):
 if 'error' in v:failures.append(dict(**item,error=v['error']));continue
 for obs,m in enumerate(v['value']['moments']):
  result=check(i['means'][obs]/i['averageQuasiDispersion'],i['dispersion'],m)
  if 'value' in old:result['priorDirectRelativeErrors']={k:abs(m[k]/old['value']['moments'][obs][k]-1) for k in ['mean','variance','devianceScale','degreesOfFreedom']}
  checks.append(dict(case=item['case'],method=item['method'],featureIndex=item['featureIndex'],observation=obs,**result))
  if result['status']!='passed' or max(result.get('priorDirectRelativeErrors',{}).values(),default=0)>2e-7:failures.append(checks[-1])
 print('checked',item['case'],item['method'],item['featureIndex'],flush=True)
(a.root/'twenty-checks.json').write_text(json.dumps(dict(checks=checks,failures=failures,status='passed' if len(checks)==124 and not failures else 'failed'),sort_keys=True,indent=2)+'\n')
grid=[dict(mean=mu,dispersion=phi) for mu in [10000,100000,1000000] for phi in [.001,.03,.1]]+[dict(mean=10000,dispersion=phi) for phi in [.7,1,4]]
payload=json.dumps(dict(moments=grid),separators=(',',':')).encode();(a.root/'extended-grid-input.json').write_bytes(payload)
assert not (a.root/'extended-grid-native.json').exists()
with (a.root/'extended-grid-native.log').open('wb') as log:
 run=subprocess.run(['ssh','macmini','/usr/bin/time','-l',a.binary],input=payload,stdout=subprocess.PIPE,stderr=log)
assert run.returncode==0;(a.root/'extended-grid-native.json').write_bytes(run.stdout)
gridchecks=[]
for i,v in zip(grid,json.loads(run.stdout)['moments'],strict=True):
 result=check(i['mean'],i['dispersion'],v['value']) if 'value' in v else dict(status='failed',error=v['error'])
 gridchecks.append(dict(input=i,**result));print('grid',i,result['status'],flush=True)
(a.root/'extended-grid-checks.json').write_text(json.dumps(dict(checks=gridchecks,status='passed' if all(r['status']=='passed' for r in gridchecks) else 'failed'),sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(twentyCheckedMoments=len(checks),twentyFailures=len(failures),gridFailures=sum(r['status']!='passed' for r in gridchecks),maximumRelativeError=max(max(r['relativeErrors'].values()) for r in checks+gridchecks if 'relativeErrors' in r),evaluations=sum(r['evaluatedCounts'] for r in checks)),indent=2))
