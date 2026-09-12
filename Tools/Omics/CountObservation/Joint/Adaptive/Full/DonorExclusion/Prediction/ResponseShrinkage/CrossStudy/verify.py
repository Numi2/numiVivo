from pathlib import Path
import json,hashlib
import numpy as np
s=Path(__file__).parent;r=json.loads((s/'prediction-freeze.json').read_text());assert r['status']=='complete' and not r['scoringStarted'];report={'folds':{}}
for tag,rec in r['folds'].items():
 for suffix,key in [('input','inputSHA256'),('output','outputSHA256')]:assert hashlib.sha256((s/'native'/(tag+'-'+suffix+'.json')).read_bytes()).hexdigest()==rec[key]
 a=json.loads((s/'native'/(tag+'-input.json')).read_text());b=json.loads((s/'native'/(tag+'-output.json')).read_text());m=b['model'];x=np.array(a['controls']).T;y=np.array(a['treated']).T;d=y-x;n=x.shape[1];assert m['featureIDs']==a['featureIDs'] and m['trainingDonorIDs']==a['trainingDonorIDs'] and a['queryDonorID'] not in m['trainingDonorIDs'];loss=[]
 def fit(xx,dd,lam):
  center=xx.mean(axis=1);mean=dd.mean(axis=1);v=xx-center[:,None];den=np.einsum('ij,ij->i',v,v)
  if lam is None:return center,mean,np.zeros_like(center)
  den+=lam*(xx.shape[1]-1);num=np.einsum('ij,ij->i',v,dd-mean[:,None]);return center,mean,np.divide(num,den,out=np.zeros_like(num),where=den>0)
 for lam in m['penalties']:
  errors=[]
  for j in range(n):
   keep=np.arange(n)!=j;c,mean,slope=fit(x[:,keep],d[:,keep],lam);p=np.maximum(0,(1+slope)*x[:,j]+mean-slope*c);errors.append(np.sum((p-y[:,j])**2))
  loss.append(float(sum(errors)/(x.shape[0]*n)))
 np.testing.assert_allclose(loss,m['innerDonorMSE'],rtol=1e-9,atol=1e-12);best=min(range(len(loss)),key=lambda i:(loss[i],-i));assert m.get('selectedPenalty')==m['penalties'][best]
 c,mean,slope=fit(x,d,m.get('selectedPenalty'))
 for key,val in [('controlMean',c),('responseMean',mean),('responseSlope',slope)]:np.testing.assert_allclose(m[key],val,rtol=1e-8,atol=1e-10)
 q=np.array(a['queryControl']);expected=np.maximum(0,(1+slope)*q+mean-slope*c);np.testing.assert_allclose(b['predictedTreated'],expected,rtol=1e-9,atol=1e-10)
 report['folds'][tag]={'features':len(q),'selectedPenalty':m.get('selectedPenalty'),'maximumPredictionDifference':float(np.max(np.abs(expected-b['predictedTreated']))),'maximumCVLossDifference':float(np.max(np.abs(np.array(loss)-m['innerDonorMSE'])))}
report['status']='pass-all-models-all-penalties-all-predictions';(s/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(report['status'],len(report['folds']))
