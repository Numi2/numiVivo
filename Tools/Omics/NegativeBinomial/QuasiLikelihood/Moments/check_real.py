#!/usr/bin/env python3
"""Independent PMF sums for every native real-data moment, in bounded batches."""
import argparse,gzip,hashlib,json
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.special import gammaln,xlogy
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ql-root',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--available',action='store_true');a=p.parse_args();sha=lambda b:hashlib.sha256(b).hexdigest();scriptHash=sha(Path(__file__).read_bytes())
manifest=json.loads((a.root/'protocol.json').read_text());results=[]

def precise(mu,phi,low,high):
 with mp.workdps(80):
  m=mp.mpf(float(mu));aa=mp.mpf(float(phi));r=1/aa;y=mp.mpf(low)
  pr=mp.exp(mp.loggamma(y+r)-mp.loggamma(r)-mp.loggamma(y+1)+y*mp.log(aa*m)-(y+r)*mp.log1p(aa*m))
  mass=[];first=[];second=[]
  for k in range(low,high+1):
   yy=mp.mpf(k);dev=2*mp.log1p(aa*m)/aa if k==0 else 2*(yy*mp.log(yy/m)-(yy+r)*mp.log((yy+r)/(m+r)))
   mass.append(pr);first.append(pr*dev);second.append(pr*dev*dev);pr*=m/(m+r)*(yy+r)/(yy+1)
  total=mp.fsum(mass);mean=mp.fsum(first)/total;var=mp.fsum(second)/total-mean*mean
  return np.asarray([float(mean),float(var),float(2*mean/var),float(2*mean*mean/var)])

for c in manifest['cases']:
 source=a.ql_root/c['id'];sourceInput=(source/'input.json.gz').read_bytes();sourceReference=(source/'reference.json.gz').read_bytes()
 original=json.loads(gzip.decompress(sourceInput));reference=json.loads(gzip.decompress(sourceReference))['results']
 for method in ['nativeTrend-adjusted','edgeRTrend-adjusted']:
  d=a.root/c['id']/method
  if not (d/'run.json').exists():
   assert a.available,(c['id'],method,'not terminal');continue
  receipt=json.loads((d/'run.json').read_text());assert sha(sourceInput)==receipt['inputSHA256'] and sha(sourceReference)==receipt['referenceSHA256']
  raw=(d/'native.json.gz').read_bytes();assert sha(raw)==receipt['outputSHA256']
  if (d/'independent-check.json').exists():
   old=json.loads((d/'independent-check.json').read_text());assert old['checkerSHA256']==scriptHash and old['outputSHA256']==receipt['outputSHA256'];results.append(old);continue
  native=json.loads(gzip.decompress(raw))['residuals'];stage=reference[method];specs=[];failures=[];refinements=[];maximum=0.;massError=0.;supportTerms=0
  for j,(v,t) in enumerate(zip(native,stage['table'],strict=True)):
   if 'error' in v:failures.append(dict(featureIndex=original['featureIndices'][j],error=v['error']));continue
   for k,mm in enumerate(v['value']['moments']):
    mu=stage['means'][j][k]/t['averageQLScale'];phi=t['actualDispersion'];assert mm['evaluatedCounts']==mm['lastCount']-mm['firstCount']+1
    specs.append((j,k,mu,phi,mm))
  geneErrors=np.zeros(len(native));start=0
  while start<len(specs):
   stop=start;work=0
   while stop<len(specs) and (work<250000 or stop==start):work+=specs[stop][4]['evaluatedCounts'];stop+=1
   batch=specs[start:stop];lengths=np.asarray([s[4]['evaluatedCounts'] for s in batch]);offsets=np.r_[0,np.cumsum(lengths)[:-1]]
   y=np.concatenate([np.arange(s[4]['firstCount'],s[4]['lastCount']+1,dtype=float) for s in batch]);mu=np.repeat([s[2] for s in batch],lengths);alpha=np.repeat([s[3] for s in batch],lengths);r=1/alpha
   logp=gammaln(y+r)-gammaln(r)-gammaln(y+1)+y*(np.log(alpha)+np.log(mu))-(y+r)*np.log1p(alpha*mu)
   logp[y==0]=-r[y==0]*np.log1p(alpha[y==0]*mu[y==0]);pr=np.exp(logp)
   dev=2*(xlogy(y,y/mu)-(y+r)*np.log1p((y-mu)/(mu+r)))
   dev[y==0]=2*np.log1p(alpha[y==0]*mu[y==0])/alpha[y==0]
   total=np.add.reduceat(pr,offsets);first=np.add.reduceat(pr*dev,offsets)/total;variance=np.add.reduceat(pr*dev*dev,offsets)/total-first*first
   expected=np.stack([first,variance,2*first/variance,2*first*first/variance],axis=1)
   supplied=np.asarray([[s[4][key] for key in ['mean','variance','devianceScale','degreesOfFreedom']] for s in batch]);errors=np.max(np.abs(supplied/expected-1),axis=1)
   for index in np.flatnonzero(~np.isfinite(errors)|(errors>2e-7)):
    j,k,m,phi,mm=batch[index];before=float(errors[index]);exact=precise(m,phi,mm['firstCount'],mm['lastCount']);expected[index]=exact;errors[index]=np.max(np.abs(supplied[index]/exact-1));refinements.append(dict(featureIndex=original['featureIndices'][j],observation=k,initialRelativeError=before,refinedRelativeError=float(errors[index])))
   for s,error,mass in zip(batch,errors,total):
    j,k,_,_,mm=s;geneErrors[j]=max(geneErrors[j],float(error));massError=max(massError,float(abs(mass-1)));maximum=max(maximum,float(error))
    if not np.isfinite(error) or error>2e-7:failures.append(dict(featureIndex=original['featureIndices'][j],observation=k,relativeMomentError=float(error)))
    if not np.isfinite(mass) or abs(mass-1)>1e-7+mm['omittedProbabilityBound']:failures.append(dict(featureIndex=original['featureIndices'][j],observation=k,independentPMFMass=float(mass)))
   supportTerms+=work;start=stop
  result=dict(case=c['id'],method=method,status='passed' if not failures else 'failed',moments=len(specs),supportTerms=supportTerms,maximumRelativeMomentError=maximum,maximumPMFMassError=massError,failures=failures,highPrecisionRefinements=refinements,checkerSHA256=scriptHash,outputSHA256=receipt['outputSHA256'])
  (d/'independent-gene-errors.json.gz').write_bytes(gzip.compress(json.dumps(dict(featureIndices=original['featureIndices'],maximumRelativeMomentErrors=geneErrors.tolist()),separators=(',',':')).encode(),mtime=0));(d/'independent-check.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');results.append(result);print(c['id'],method,result['status'],len(specs),maximum,flush=True)
summary=dict(status='passed' if len(results)==58 and all(r['status']=='passed' for r in results) else 'incomplete-or-failed',arms=len(results),moments=sum(r['moments'] for r in results),supportTerms=sum(r['supportTerms'] for r in results),maximumRelativeMomentError=max((r['maximumRelativeMomentError'] for r in results),default=0),maximumPMFMassError=max((r['maximumPMFMassError'] for r in results),default=0),refinements=sum(len(r['highPrecisionRefinements']) for r in results),failures=[r for r in results if r['status']!='passed'],checkerSHA256=scriptHash)
(a.root/('independent-available-summary.json' if a.available else 'independent-summary.json')).write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='failures'},indent=2))
