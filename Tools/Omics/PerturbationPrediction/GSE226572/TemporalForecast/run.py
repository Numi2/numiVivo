from pathlib import Path
import numpy as np,json,hashlib,subprocess,time,math
r=Path(__file__).parent;source=Path('/Users/n/numivivo-duration-source-20260912.npz');h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();assert h(source)=='534a0f5bcd676af37cda6af9b7085da012dc68edc90ac99c3fdffc06b971e341'
protocol={'scope':'Reused GSE226572 development endpoint; full QC-admitted PBMC population, not B-cell prediction','split':'Exclude entire query donor and all times greater than or equal to query time; all 18 folds retained','model':'Forecast from equal-donor mean response at latest earlier training time; no eligible treated observation falls back to no change','baseline':'equal-donor, equal-within-donor-time mean response from identical training subset; plus no change','features':12993,'gate':'at least 5 percent lower equal-within-donor-time mean RMSE against both baselines in each of three held-out donors','sourceSHA256':h(source),'nativeSourceSHA256':h(r/'Temporal.swift'),'createdUnix':time.time()}
(r/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');a=np.load(source);donors=['D34','D38','D39'];obs=[]
for donor in donors:
 prefix='query:GSE226572:'+donor;control=a[prefix+'|noChange']
 for key in a.files:
  if key.startswith(prefix+'|truth|IFNB:'):obs.append((donor,float(key.rsplit(':',1)[1][:-1]),a[key],control))
assert len(obs)==18 and len(a['featureIDs'])==12993
(r/'inputs').mkdir();(r/'outputs').mkdir();records=[]
for i,(donor,t,target,control) in enumerate(obs):
 training=[{'donor':d,'time':v,'response':(y-c).tolist()} for d,v,y,c in obs if d!=donor and v<t]
 data={'donor':donor,'time':t,'control':control.tolist(),'training':training};inp=r/'inputs'/f'{i}.json';out=r/'outputs'/f'{i}.json';inp.write_text(json.dumps(data,separators=(',',':')))
 start=time.monotonic();subprocess.run([str(r/'temporal'),str(inp),str(out)],check=True);records.append({'fold':i,'donor':donor,'time':t,'inputSHA256':h(inp),'outputSHA256':h(out),'seconds':time.monotonic()-start})
(r/'freeze.json').write_text(json.dumps({'scored':False,'binarySHA256':h(r/'temporal'),'records':records,'createdUnix':time.time()},indent=2)+'\n')
rows=[];maxerror=0
for i,(donor,t,target,control) in enumerate(obs):
 b=json.loads((r/'outputs'/f'{i}.json').read_text());training=[(d,v,y-c) for d,v,y,c in obs if d!=donor and v<t];ts=sorted({v for d,v,y in training});knots=[0]+ts;curve=np.vstack([np.zeros_like(control)]+[np.mean([y for d,v,y in training if v==time],axis=0) for time in ts]);prediction=np.maximum(0,control+np.array([np.interp(t,knots,curve[:,j]) for j in range(len(control))]));baseline=np.maximum(0,control+np.mean([np.mean([y for d,v,y in training if d==who],axis=0) for who in sorted({d for d,v,y in training})],axis=0)) if training else control.copy()
 for ref,key in [(prediction,'predicted'),(baseline,'baseline')]:
  err=float(np.max(np.abs(ref-b[key])));maxerror=max(maxerror,err);assert err<1e-12
 scores={}
 for key,v in [('candidate',b['predicted']),('trainingMean',b['baseline']),('noChange',control)]:
  value=float(np.sqrt(np.mean((np.array(v)-target)**2)));scalar=math.sqrt(math.fsum((float(x)-float(y))**2 for x,y in zip(v,target))/len(target));assert abs(value-scalar)<1e-12;scores[key]=value
 rows.append({'donor':donor,'time':t,'trainingTimes':ts,'noEarlierTreatedObservation':not bool(training),'lowerTime':b['lowerTime'],'upperTime':b['upperTime'],'rmse':scores})
summary={}
for donor in donors:
 chosen=[x for x in rows if x['donor']==donor];means={k:float(np.mean([x['rmse'][k] for x in chosen])) for k in ['candidate','trainingMean','noChange']};gains={k:1-means['candidate']/means[k] for k in ['trainingMean','noChange']};summary[donor]={'meanRMSE':means,'gainFraction':gains,'pass':all(x>=.05 for x in gains.values())}
(r/'results.json').write_text(json.dumps({'rows':rows,'donors':summary,'maximumNativeReferenceDifference':maxerror,'allDonorsPass':all(x['pass'] for x in summary.values()),'predictedValues':18*12993},indent=2)+'\n');print(json.dumps(summary))
