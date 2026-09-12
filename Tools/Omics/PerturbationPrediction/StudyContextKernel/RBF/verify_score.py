from pathlib import Path
import json,hashlib,math
import numpy as np
s=Path(__file__).parent;read=lambda p:json.loads(p.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
e=read(s/'execution.json');f=read(s/'prediction-freeze.json');assert sha(s/'execution.json')==f['executionSHA256'];assert sha(s/'ContextKernel.swift')==e['sourceSHA256'];assert sha(s/'protocol.json')==e['protocolSHA256']
source=s.parent/'numivivo-parse-context-evaluation-20260912';data=np.load(source/'source-cohorts.npz')
scales=[1];penalties=[.01,.1,1,10,100,None];grid=[(z,l) for z in scales for l in penalties]
def ref(x,y,q,st,l,z):
 w=np.array([1/(len(set(st))*st.count(v)) for v in st]);g=x.shape[1]
 d=np.maximum(0,(np.sum(x*x,axis=1)[:,None]+np.sum(x*x,axis=1)[None,:]-2*x@x.T)/g)
 bandwidth=float(np.median(d[np.tril_indices(len(x),-1)]));bandwidth=bandwidth if bandwidth>0 else 1.0
 raw=np.exp(-d/(2*bandwidth));means=raw@w;grand=w@means
 if l is None:beta=np.broadcast_to(w,(len(q),len(w))).copy()
 else:
  k=(raw-means[:,None]-means[None,:]+grand)*np.sqrt(w[:,None]*w[None,:])
  dq=np.maximum(0,(np.sum(x*x,axis=1)[:,None]+np.sum(q*q,axis=1)[None,:]-2*x@q.T)/g)
  rq=np.exp(-dq/(2*bandwidth));rhs=(rq-means[:,None]-(w@rq)[None,:]+grand)*np.sqrt(w[:,None])
  ev,u=np.linalg.eigh(k);a=u@((u.T@rhs)/(ev[:,None]+l));b=a.T*np.sqrt(w);beta=w+b-b.sum(1)[:,None]*w

 return np.maximum(0,q+z*(beta@(y-x))),beta
out={}
for name,h in e['inputs'].items():
 assert sha(s/'inputs'/name)==h and sha(s/'native'/name)==f['outputs'][name]['SHA256']
 a=read(s/'inputs'/name);b=read(s/'native'/name);x=np.array(a['controls']);y=np.array(a['treated']);q=np.array(a['queryControls']);st=a['trainingStudies'];loss=[]
 assert a['featureIDs']==b['featureIDs'] and a['queryDonorIDs']==b['queryDonorIDs']
 for study in sorted(set(st)):
  keep=np.array([v!=study for v in st]);test=~keep;row=[]
  for z,l in grid:
   pred,_=ref(x[keep],y[keep],x[test],[v for v,k in zip(st,keep) if k],l,z);row.append(float(np.mean((pred-y[test])**2)))
  loss.append(row)
 means=np.mean(loss,axis=0);best=min(range(len(grid)),key=lambda i:(means[i],-i));z,l=grid[best]
 assert b['model'].get('selectedPenalty')==l
 np.testing.assert_allclose(loss,b['model']['innerStudyMSE'],rtol=1e-10,atol=1e-11)
 pred,beta=ref(x,y,q,st,l,z);baseline,_=ref(x,y,q,st,None,1)
 np.testing.assert_allclose(pred,b['predictedTreated'],rtol=1e-10,atol=1e-10);np.testing.assert_allclose(beta,b['model']['queryResponseWeights'],rtol=1e-9,atol=1e-10);np.testing.assert_allclose(baseline,b['trainingMeanPredictions'],rtol=1e-12,atol=1e-12)
 target=data[a['queryStudy']+'_treated'];rows=[];maxscore=0
 for i,d in enumerate(a['queryDonorIDs']):
  scores={}
  for key,v in [('candidate',b['predictedTreated'][i]),('trainingMean',b['trainingMeanPredictions'][i]),('noChange',q[i])]:
   rmse=float(np.sqrt(np.mean((np.array(v)-target[i])**2)));scalar=math.sqrt(math.fsum((float(k)-float(t))**2 for k,t in zip(v,target[i]))/len(v));assert abs(rmse-scalar)<1e-12;maxscore=max(maxscore,abs(rmse-scalar));scores[key]=rmse
  rows.append(dict(donor=d,rmse=scores))
 avg={k:math.fsum(r['rmse'][k] for r in rows)/len(rows) for k in scores};gains={k:1-avg['candidate']/avg[k] for k in ['trainingMean','noChange']}
 out[a['queryStudy']]=dict(selectedScale=z,selectedPenalty=l,donors=rows,meanRMSE=avg,gainFraction=gains,status='PASS' if all(v>=.05 for v in gains.values()) else 'FAIL',maximumPredictionDifference=float(np.max(abs(pred-b['predictedTreated']))),maximumScoreDifference=maxscore)
 print(a['queryStudy'],z,l,gains,out[a['queryStudy']]['status'],flush=True)
result=dict(scope='Reused development studies only. No Parse evaluation.',studies=out,protocolSHA256=sha(s/'protocol.json'),predictionFreezeSHA256=sha(s/'prediction-freeze.json'),verifierSHA256=sha(Path(__file__)),status='PASS' if all(v['status']=='PASS' for v in out.values()) else 'FAIL')
with (s/'verified-scores.json').open('x') as f:json.dump(result,f,indent=2)
