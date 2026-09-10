#!/usr/bin/env python3
"""Check native moment/tail results against independently summed probability masses."""
import argparse,gzip,hashlib,json,subprocess,time
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.stats import nbinom,poisson
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--binary',required=True);p.add_argument('--host',default='macmini');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
sha=lambda b:hashlib.sha256(b).hexdigest();cases=[dict(mean=mu,dispersion=phi) for mu in [1e-8,1e-4,.01,.1,.5,1,3,10,30,100,1000] for phi in [0,1e-8,.001,.1,.7,1,4]]
assert not (a.out/'grid-native.json.gz').exists();payload=json.dumps(dict(moments=cases),separators=(',',':')).encode();(a.out/'grid-input.json').write_bytes(payload)
t0=time.monotonic();run=subprocess.run(['ssh',a.host,a.binary],input=payload,stdout=subprocess.PIPE,stderr=subprocess.PIPE);elapsed=time.monotonic()-t0
(a.out/'grid-native.log').write_bytes(run.stderr);assert run.returncode==0;packed=gzip.compress(run.stdout,mtime=0);(a.out/'grid-native.json.gz').write_bytes(packed);out=json.loads(run.stdout)['moments'];rows=[]
for i,(c,v) in enumerate(zip(cases,out,strict=True)):
 mu,phi=c['mean'],c['dispersion'];dist=poisson(mu) if phi==0 else nbinom(1/phi,1/(1+phi*mu));upper=int(dist.ppf(1-1e-15))+128;k=np.arange(upper+1,dtype=float)
 pr=dist.pmf(k);np.testing.assert_allclose(pr.sum(),1,atol=1e-8,rtol=0)
 if mu<=.01 or phi<=1e-8:
  with mp.workdps(80):
   m=mp.mpf(float(mu));alpha=mp.mpf(float(phi));r=1/alpha if phi else None
   prob=mp.exp(-m) if phi==0 else mp.exp(-r*mp.log1p(alpha*m));terms=[];terms2=[]
   for count in range(upper+1):
    y=mp.mpf(count)
    d=2*m if count==0 and phi==0 else 2*mp.log1p(alpha*m)/alpha if count==0 else 2*(y*mp.log(y/m)-y+m) if phi==0 else 2*(y*mp.log(y/m)-(y+r)*mp.log((y+r)/(m+r)))
    terms.append(prob*d);terms2.append(prob*d*d)
    prob*=m/(count+1) if phi==0 else (m/(m+r))*(count+r)/(count+1)
   mean=mp.fsum(terms);var=mp.fsum(terms2)-mean*mean;expected=dict(mean=float(mean),variance=float(var),devianceScale=float(2*mean/var),degreesOfFreedom=float(2*mean*mean/var))
  owner='80-digit direct probability recurrence'
 else:
  with np.errstate(divide='ignore',invalid='ignore'):
   first=np.where(k==0,0,k*np.log(k/mu));d=2*(first-(k+1/phi)*np.log((k+1/phi)/(mu+1/phi)))
  mean=float(pr@d);var=float(pr@(d*d)-mean*mean);expected=dict(mean=mean,variance=var,devianceScale=2*mean/var,degreesOfFreedom=2*mean*mean/var);owner='SciPy PMF and direct saturated likelihood'
 errors=None if 'value' not in v else {key:abs(v['value'][key]/value-1) for key,value in expected.items()}
 row=dict(index=i,input=c,reference=expected,referenceOwner=owner,referenceUpperCount=upper,native=v,relativeErrors=errors,status='passed' if errors and max(errors.values())<=2e-7 else 'failed');rows.append(row)
result=dict(cases=len(rows),status='passed' if all(r['status']=='passed' for r in rows) else 'failed',seconds=elapsed,maximumRelativeError=max(max(r['relativeErrors'].values()) for r in rows if r['relativeErrors']),maximumEvaluatedCounts=max(r['native']['value']['evaluatedCounts'] for r in rows if 'value' in r['native']),binarySHA256=subprocess.check_output(['ssh',a.host,'shasum','-a','256',a.binary]).decode().split()[0],inputSHA256=sha(payload),outputSHA256=sha(packed),protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md').read_bytes()),rows=rows)
(a.out/'grid-check.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
