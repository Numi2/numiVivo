#!/usr/bin/env python3
"""Run the predeclared higher-work sensitivity and independently sum every moment."""
import argparse,gzip,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from scipy.special import gammaln,xlogy
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ql-root',type=Path,required=True);p.add_argument('--initial-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--binary',required=True);a=p.parse_args();sha=lambda b:hashlib.sha256(b).hexdigest()
manifest=json.loads((a.initial_root/'protocol.json').read_text());runs=json.loads((a.initial_root/'complete.json').read_text());items=[];rows=[]
for run in runs:
 if not run['failures']:continue
 c=run['case'];source=a.ql_root/c;raw=(source/'input.json.gz').read_bytes();assert sha(raw)==run['inputSHA256'];i=json.loads(gzip.decompress(raw));raw=(source/'reference.json.gz').read_bytes();assert sha(raw)==run['referenceSHA256'];stage=json.loads(gzip.decompress(raw))['results'][run['method']]
 for fail in run['failures']:
  assert 'exhausted 1000000' in fail['error'];j=i['featureIndices'].index(fail['featureIndex']);t=stage['table'][j]
  rows.append(dict(counts=i['counts'][j],means=stage['means'][j],design=i['design'],dispersion=t['actualDispersion'],averageQuasiDispersion=t['averageQLScale']))
  items.append(dict(case=c,method=run['method'],featureIndex=fail['featureIndex'],originalError=fail['error'],inputSHA256=run['inputSHA256'],referenceSHA256=run['referenceSHA256']))
assert len(rows)==20;a.out.mkdir(parents=True,exist_ok=True);assert not (a.out/'native.json.gz').exists();raw=json.dumps(rows,separators=(',',':')).encode();(a.out/'input.json').write_bytes(raw)
metadata=dict(items=items,maximumTerms=10000000,relativeTolerance=1e-10,protocolSHA256=sha(Path(__file__).with_name('WORK_LIMIT_SENSITIVITY.md').read_bytes()),requestSHA256=sha(raw),binarySHA256=subprocess.check_output(['ssh',manifest['host'],'shasum','-a','256',a.binary]).decode().split()[0]);(a.out/'protocol.json').write_text(json.dumps(metadata,sort_keys=True,indent=2)+'\n')
t0=time.monotonic()
with (a.out/'native.log').open('wb') as log:r=subprocess.run(['ssh',manifest['host'],'/usr/bin/time','-l',a.binary],input=raw,stdout=subprocess.PIPE,stderr=log)
assert r.returncode==0;packed=gzip.compress(r.stdout,mtime=0);(a.out/'native.json.gz').write_bytes(packed);native=json.loads(r.stdout);elapsed=time.monotonic()-t0;result=[];failures=[]
for item,input_,v in zip(items,rows,native,strict=True):
 if 'error' in v:failures.append(dict(**item,error=v['error']));continue
 for obs,m in enumerate(v['value']['moments']):
  mu=input_['means'][obs]/input_['averageQuasiDispersion'];alpha=input_['dispersion'];size=1/alpha;mass=first=second=0.
  for start in range(m['firstCount'],m['lastCount']+1,250000):
   k=np.arange(start,min(start+250000,m['lastCount']+1),dtype=float);pr=np.exp(gammaln(k+size)-gammaln(size)-gammaln(k+1)+k*(np.log(alpha)+np.log(mu))-(k+size)*np.log1p(alpha*mu));pr[k==0]=np.exp(-size*np.log1p(alpha*mu))
   d=2*(xlogy(k,k/mu)-(k+size)*np.log1p((k-mu)/(mu+size)));d[k==0]=2*np.log1p(alpha*mu)/alpha
   mass+=float(pr.sum());first+=float(pr@d);second+=float(pr@(d*d))
  mean=first/mass;var=second/mass-mean*mean;expected=dict(mean=mean,variance=var,devianceScale=2*mean/var,degreesOfFreedom=2*mean*mean/var);errors={k:abs(m[k]/x-1) for k,x in expected.items()}
  entry=dict(**item,observation=obs,evaluatedCounts=m['evaluatedCounts'],relativeErrors=errors,referenceMass=mass,reference=expected)
  result.append(entry)
  if max(errors.values())>2e-7 or abs(mass-1)>1e-7+m['omittedProbabilityBound']:failures.append(entry)
summary=dict(status='passed' if not failures and len(native)==20 else 'failed',genes=len(native),moments=len(result),evaluatedCounts=sum(v['evaluatedCounts'] for v in result),maximumEvaluatedCounts=max((v['evaluatedCounts'] for v in result),default=0),maximumRelativeMomentError=max((max(v['relativeErrors'].values()) for v in result),default=0),maximumPMFMassError=max((abs(v['referenceMass']-1) for v in result),default=0),seconds=elapsed,outputSHA256=sha(packed),logicalSHA256=sha(r.stdout),binarySHA256=metadata['binarySHA256'],failures=failures,checks=result)
(a.out/'summary.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k not in ['checks','failures']},indent=2))
