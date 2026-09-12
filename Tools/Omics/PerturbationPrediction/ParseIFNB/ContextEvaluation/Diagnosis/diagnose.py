from pathlib import Path
import json,hashlib,math
import numpy as np
s=Path(__file__).parent;o=s.parent/'numivivo-parse-context-diagnosis-20260912';o.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
f=read(s/'prediction-freeze.json');scores=read(s/'scores.json');assert sha(s/'native/Parse.json')==f['outputSHA256'] and sha(s/'targets.npy')==scores['targetsSHA256']
a=read(s/'inputs/Parse.json');b=read(s/'native/Parse.json');t=np.load(s/'targets.npy');q=np.array(a['queryControls']);observed=t-q
rows=[]
for name,key in [('candidate','predictedTreated'),('trainingMean','trainingMeanPredictions')]:
 response=np.array(b[key])-q
 for i,donor in enumerate(a['queryDonorIDs']):
  r=response[i];y=observed[i];rr=float(r@r);yy=float(y@y);ry=float(r@y);g=len(r)
  # Error decomposition uses applied response after native nonnegative clipping.
  error=(rr+yy-2*ry)/g;direct=float(np.mean((r-y)**2));assert abs(error-direct)<1e-12
  scalar=math.fsum((float(x)-float(z))**2 for x,z in zip(r,y))/g;assert abs(scalar-direct)<1e-12
  rows.append(dict(model=name,donor=donor,predictedResponseRMS=math.sqrt(rr/g),observedResponseRMS=math.sqrt(yy/g),responseNormRatio=math.sqrt(rr/yy) if yy else None,cosine=ry/math.sqrt(rr*yy) if rr*yy else None,responsePower=rr/g,alignmentCredit=2*ry/g,excessMSEOverNoChange=(rr-2*ry)/g,rmse=math.sqrt(direct)))
out=dict(scope='Post-hoc diagnosis only. No model refit, scale optimization, donor exclusion or new validation claim.',definition='Applied predicted response = clipped predicted treated minus query control. Excess MSE over no change = response power minus twice observed/predicted inner product per feature.',sources={str(s/n):sha(s/n) for n in ['inputs/Parse.json','native/Parse.json','targets.npy','scores.json','prediction-freeze.json']},scriptSHA256=sha(Path(__file__)),rows=rows)
(o/'diagnosis.json').write_text(json.dumps(out,indent=2)+'\n')
for name in ['candidate','trainingMean']:
 v=[x for x in rows if x['model']==name];print(name,dict(positiveCosine=sum(x['cosine']>0 for x in v),normRatioRange=[min(x['responseNormRatio'] for x in v),max(x['responseNormRatio'] for x in v)],cosineRange=[min(x['cosine'] for x in v),max(x['cosine'] for x in v)],worseThanNoChange=sum(x['excessMSEOverNoChange']>0 for x in v)))
