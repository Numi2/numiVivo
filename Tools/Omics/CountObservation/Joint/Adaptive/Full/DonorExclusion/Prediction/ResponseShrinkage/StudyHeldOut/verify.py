from pathlib import Path
import json,hashlib
import numpy as np
s=Path(__file__).parent;fr=json.loads((s/'prediction-freeze.json').read_text());assert fr['status']=='complete-three-study-native-predictions' and not fr['scoringStarted'];out={'folds':{}}
def weights(studies):
 unique=set(studies);return np.array([1/(len(unique)*studies.count(v)) for v in studies])
def fit(x,d,w,lam):
 c=x@w;m=d@w;xc=x-c[:,None];den=np.sum(xc**2*w,axis=1)
 if lam is None:return c,m,np.zeros_like(c)
 den+=lam*(1-np.sum(w*w));cov=np.sum(xc*(d-m[:,None])*w,axis=1);return c,m,np.divide(cov,den,out=np.zeros_like(cov),where=den>0)
for study,rec in fr['folds'].items():
 for folder,key in [('inputs','inputSHA256'),('native','outputSHA256')]:assert hashlib.sha256((s/folder/(study+'.json')).read_bytes()).hexdigest()==rec[key]
 a=json.loads((s/'inputs'/(study+'.json')).read_text());b=json.loads((s/'native'/(study+'.json')).read_text());m=b['model'];ss=a['trainingStudies'];unique=sorted(set(ss));assert study not in unique and b['queryStudy']==study and m['trainingStudies']==unique and m['trainingDonorIDs']==a['trainingDonorIDs'] and b['queryDonorIDs']==a['queryDonorIDs'] and set(a['trainingDonorIDs']).isdisjoint(a['queryDonorIDs']) and b['featureIDs']==a['featureIDs'];x=np.array(a['controls']).T;y=np.array(a['treated']).T;d=y-x;loss=[]
 for held in unique:
  train=np.array([i for i,v in enumerate(ss) if v!=held]);test=np.array([i for i,v in enumerate(ss) if v==held]);w=weights([ss[i] for i in train]);errors=[]
  for lam in m['penalties']:
   c,mean,slope=fit(x[:,train],d[:,train],w,lam);p=np.maximum(0,(1+slope[:,None])*x[:,test]+(mean-slope*c)[:,None]);errors.append(float(np.mean((p-y[:,test])**2)))
  loss.append(errors)
 np.testing.assert_allclose(loss,m['innerStudyMSE'],rtol=1e-9,atol=1e-11);means=np.mean(loss,axis=0);best=min(range(len(means)),key=lambda i:(means[i],-i));assert m.get('selectedPenalty')==m['penalties'][best]
 w=weights(ss);np.testing.assert_allclose(w,m['donorWeights'],rtol=1e-12,atol=1e-12)
 for group in unique:assert abs(w[np.array(ss)==group].sum()-1/len(unique))<1e-12
 c,mean,slope=fit(x,d,w,m.get('selectedPenalty'))
 for key,val in [('controlMean',c),('responseMean',mean),('responseSlope',slope)]:np.testing.assert_allclose(m[key],val,rtol=1e-8,atol=1e-10)
 q=np.array(a['queryControls']);p=np.maximum(0,q*(1+slope)+mean-slope*c);base=np.maximum(0,q+mean);np.testing.assert_allclose(p,b['predictedTreated'],rtol=1e-9,atol=1e-10);np.testing.assert_allclose(base,b['trainingMeanPredictions'],rtol=1e-9,atol=1e-10)
 out['folds'][study]={'queryDonors':len(q),'features':len(c),'selectedPenalty':m.get('selectedPenalty'),'maximumPredictionDifference':float(np.max(np.abs(p-b['predictedTreated']))),'maximumInnerLossDifference':float(np.max(np.abs(np.array(loss)-m['innerStudyMSE']))),'eachTrainingStudyWeight':1/len(unique)};print(study,out['folds'][study],flush=True)
out['status']='pass-all-weights-models-inner-study-losses-and-predictions';(s/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
