#!/usr/bin/env python3
"""Compare the complete adaptive family with retained direct and independent evidence."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--old',type=Path,required=True);p.add_argument('--ql-root',type=Path,required=True);a=p.parse_args();sha=lambda b:hashlib.sha256(b).hexdigest()
newRuns=json.loads((a.root/'family/complete.json').read_text());assert len(newRuns)==58
oldRuns={(r['case'],r['method']):r for r in json.loads((a.old/'final/complete.json').read_text())}
assert {(r['case'],r['method']) for r in newRuns}==set(oldRuns) and len(oldRuns)==58
items=json.loads((a.old/'work-sensitivity/protocol.json').read_text())['items'];supp=json.loads((a.root/'twenty-native.json').read_text())['residuals']
supplement={(i['case'],i['method'],i['featureIndex']):o for i,o in zip(items,supp,strict=True)}
assert json.loads((a.root/'twenty-checks.json').read_text())['status']=='passed'
assert json.loads((a.root/'extended-grid-checks.json').read_text())['status']=='passed'
fields=['mean','variance','devianceScale','degreesOfFreedom'];results=[];failures=[]
for run in newRuns:
 c,method=run['case'],run['method'];oldRun=oldRuns[(c,method)]
 for key in ['inputSHA256','referenceSHA256','requestSHA256']:assert run[key]==oldRun[key],(c,method,key)
 raw=(a.root/'family'/c/method/'native.json.gz').read_bytes();assert sha(raw)==run['outputSHA256'];new=json.loads(gzip.decompress(raw))['residuals']
 raw=(a.old/'final'/c/method/'native.json.gz').read_bytes();assert sha(raw)==oldRun['outputSHA256'];old=json.loads(gzip.decompress(raw))['residuals']
 raw=(a.ql_root/c/'input.json.gz').read_bytes();assert sha(raw)==run['inputSHA256'];features=json.loads(gzip.decompress(raw))['featureIndices']
 errors={k:0. for k in fields};direct=adaptive=added=identical=0;observations=0;geneChecks=[]
 for feature,nov,prev in zip(features,new,old,strict=True):
  if 'error' in nov:failures.append(dict(case=c,method=method,featureIndex=feature,error=nov['error']));continue
  if 'error' in prev:
   assert nov==supplement[(c,method,feature)],(c,method,feature,'supplement differs')
   added+=len(nov['value']['moments']);prev=supplement[(c,method,feature)]
  perGene={k:0. for k in fields};nv,pv=nov['value'],prev['value']
  assert nv['leverage']==pv['leverage']
  for n,o in zip(nv['moments'],pv['moments'],strict=True):
   observations+=1;adaptive+=int('summation' in n);direct+=int('summation' not in n);identical+=int(n==o)
   for key in fields:
    relative=abs(n[key]/o[key]-1);errors[key]=max(errors[key],relative);perGene[key]=max(perGene[key],relative)
   d=n.get('summation',{});bound=max((n['meanTruncationBound']+d.get('meanErrorBound',0))/n['mean'],(n['varianceTruncationBound']+d.get('varianceErrorBound',0))/n['variance'])
   assert bound<=1e-10,(c,method,feature,bound)
  if max(perGene.values())>2e-7:failures.append(dict(case=c,method=method,featureIndex=feature,relativeErrors=perGene))
  residualErrors={k:abs(nv[k]-pv[k])/max(abs(pv[k]),1e-30) for k in ['deviance','degreesOfFreedom']}
  if max(residualErrors.values())>2e-7:failures.append(dict(case=c,method=method,featureIndex=feature,residualErrors=residualErrors))
  geneChecks.append(dict(featureIndex=feature,relativeMomentErrors=perGene,residualRelativeErrors=residualErrors))
 dest=a.root/'family'/c/method/'comparison.json.gz';dest.write_bytes(gzip.compress(json.dumps(geneChecks,separators=(',',':')).encode(),mtime=0))
 results.append(dict(case=c,method=method,moments=observations,directMoments=direct,adaptiveMoments=adaptive,newlyAvailableMoments=added,identicalMomentRecords=identical,maximumRelativeErrors=errors,comparisonSHA256=sha(dest.read_bytes()),nativeOutputSHA256=run['outputSHA256'],priorOutputSHA256=oldRun['outputSHA256']))
 print('compared',c,method,'moments',observations,'adaptive',adaptive,flush=True)
summary=dict(status='passed' if not failures and len(results)==58 else 'failed',arms=len(results),moments=sum(r['moments'] for r in results),adaptiveMoments=sum(r['adaptiveMoments'] for r in results),newlyAvailableMoments=sum(r['newlyAvailableMoments'] for r in results),maximumRelativeError=max(max(r['maximumRelativeErrors'].values()) for r in results),evaluatedCounts=sum(r['evaluatedCounts'] for r in newRuns),maximumRelativeCombinedBound=max(r['maximumRelativeCombinedBound'] for r in newRuns),failures=failures,comparisons=results)
(a.root/'family-comparison.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k not in ['comparisons','failures']},indent=2))
