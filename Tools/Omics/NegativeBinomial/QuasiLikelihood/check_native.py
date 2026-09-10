#!/usr/bin/env python3
"""Measure native unit deviances on all reference fits and exact boundary cases."""
import argparse,gzip,hashlib,json,subprocess
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.special import xlogy
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--binary',required=True);a=p.parse_args();root=a.root
sha=lambda b:hashlib.sha256(b).hexdigest()
protocol=json.loads((root/'protocol.json').read_text());host=protocol['host'];binaryHash=subprocess.check_output(['ssh',host,'shasum','-a','256',a.binary]).decode().split()[0]
def execute(rows,path):
 raw=json.dumps(rows,allow_nan=False,separators=(',',':')).encode();process=subprocess.run(['ssh',host,a.binary],input=raw,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 (path.with_suffix('.log')).write_bytes(process.stderr);assert process.returncode==0
 packed=gzip.compress(process.stdout,mtime=0);path.write_bytes(packed)
 return json.loads(process.stdout),dict(inputSHA256=sha(raw),outputSHA256=sha(packed),logicalSHA256=sha(process.stdout),binarySHA256=binaryHash)
def precise(y,mu,alpha):
 with mp.workdps(100):
  yy=mp.mpf(int(y));mm=mp.mpf(float(mu));aa=mp.mpf(float(alpha))
  if yy==mm:return 0.0
  if yy==0:return float(2*mm if aa==0 else 2*mp.log1p(aa*mm)/aa)
  return float(2*(yy*mp.log(yy/mm)-yy+mm) if aa==0 else 2*(yy*mp.log(yy/mm)-(yy+1/aa)*mp.log((yy+1/aa)/(mm+1/aa))))
# Fixed Cartesian domain probe, independent of any measured errors.
probes=[]
for y in [0,1,3,100,10**6,10**9,2**53]:
 for mu in [1e-200,1e-8,.1,1.,4.,100.,1e6,1e12,1e200]+([float(y)*v for v in [.799999999999,.8,.800000000001,.9999999999,1.,1.0000000001,1.333333333332,4/3,1.333333333334]] if y else []):
  for alpha in [0.,1e-8,1e-4,.1,1.,100.]:probes.append(dict(counts=[y],means=[mu],dispersion=alpha))
assert not (root/'boundary-native.json.gz').exists()
values,receipt=execute(probes,root/'boundary-native.json.gz');boundary=[];maxRelative=0.
for i,(row,result) in enumerate(zip(probes,values,strict=True)):
 expected=precise(row['counts'][0],row['means'][0],row['dispersion']);got=result.get('values',[None])[0]
 error=None if got is None else abs(got-expected);relative=None if got is None else error/max(abs(expected),1e-300)
 if relative is not None:maxRelative=max(maxRelative,relative)
 # Relative precision for nonzero boundary deviances; absolute error cannot
 # conceal cancellation of a tiny but representable deviance.
 if got is None or relative>2e-10:boundary.append(dict(index=i,input=row,result=result,expected=expected,relativeError=relative))
(root/'boundary-check.json').write_text(json.dumps(dict(**receipt,cases=len(probes),maximumRelativeError=maxRelative,failures=boundary),sort_keys=True,indent=2)+'\n')
results=[]
for c in protocol['cases']:
 d=root/c['id'];i=json.loads(gzip.decompress((d/'input.json.gz').read_bytes()));r=json.loads(gzip.decompress((d/'reference.json.gz').read_bytes()));rows=[];slices=[]
 for name,stage in r['results'].items():
  start=len(rows)
  rows.extend(dict(counts=y,means=mu,dispersion=t['fittedDispersion']) for y,mu,t in zip(i['counts'],stage['means'],stage['table'],strict=True));slices.append((name,start,len(rows)))
 assert not (d/'native.json.gz').exists();native,receipt=execute(rows,d/'native.json.gz');details=[];failures=[];refinements=[];maximum=0.;rawMaximum=0.;unitCount=0
 for name,start,end in slices:
  selected=rows[start:end];stage=r['results'][name];y=np.asarray([v['counts'] for v in selected],dtype=np.longdouble);mu=np.asarray([v['means'] for v in selected],dtype=np.longdouble);aa=np.asarray([v['dispersion'] for v in selected],dtype=np.longdouble)[:,None]
  # Independent direct saturated likelihood subtraction; cancellation-prone
  # comparisons are resolved at exact saved means with 100 digits.
  with np.errstate(divide='ignore',invalid='ignore'):
   logRatio=np.where(y==0,0,np.log(y/mu));reference=2*(y*logRatio-(y+1/aa)*np.log((y+1/aa)/(mu+1/aa)))
  nativeRows=native[start:end]
  for j,v in enumerate(nativeRows):
   if v.get('error'):failures.append(dict(method=name,featureIndex=i['featureIndices'][j],error=v['error']))
  if any(v.get('error') for v in nativeRows):continue
  got=np.asarray([v['values'] for v in nativeRows]);reference=np.asarray(reference,dtype=float);difference=np.abs(got-reference)
  for j,k in zip(*np.where(difference>1e-8)):
   previous=float(reference[j,k]);reference[j,k]=precise(y[j,k],mu[j,k],aa[j,0]);refinements.append(dict(method=name,featureIndex=i['featureIndices'][j],observation=int(k),initialReference=previous,preciseReference=float(reference[j,k]),native=float(got[j,k])))
  difference=np.abs(got-reference);unitCount+=got.size;maximum=max(maximum,float(difference.max()))
  bad=~np.isfinite(got)|(difference>2e-7+2e-10*np.abs(reference))
  for j,k in zip(*np.where(bad)):failures.append(dict(method=name,featureIndex=i['featureIndices'][j],observation=int(k),native=float(got[j,k]),reference=float(reference[j,k])))
  rawError=np.abs(got.sum(axis=1)-np.asarray([t['rawDeviance'] for t in stage['table']]))
  rawMaximum=max(rawMaximum,float(rawError.max()))
  details.append(dict(method=name,units=int(got.size),maximumAbsoluteError=float(difference.max()),maximumEdgeRRawDevianceDifference=float(rawError.max())))
 result=dict(case=c['id'],status='passed' if not failures else 'failed',**receipt,units=unitCount,maximumAbsoluteError=maximum,maximumEdgeRRawDevianceDifference=rawMaximum,failures=failures,precisionRefinements=refinements,stages=details)
 (d/'native-check.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');results.append(result);print(c['id'],result['status'],unitCount,maximum,flush=True)
summary=dict(status='passed' if not boundary and all(r['status']=='passed' for r in results) else 'failed',binarySHA256=binaryHash,boundaryCases=len(probes),boundaryFailures=boundary,boundaryMaximumRelativeError=maxRelative,cases=len(results),units=sum(r['units'] for r in results),maximumAbsoluteError=max(r['maximumAbsoluteError'] for r in results),maximumEdgeRRawDevianceDifference=max(r['maximumEdgeRRawDevianceDifference'] for r in results),refinements=sum(len(r['precisionRefinements']) for r in results),failures=[r for r in results if r['failures']])
(root/'native-summary.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k not in ['boundaryFailures','failures']},indent=2))
