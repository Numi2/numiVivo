from pathlib import Path
import json,hashlib,subprocess,math
import numpy as np,anndata
r=Path(__file__).parent;h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();source=Path('/Users/n/numivivo-duration-source-20260912.npz');assert h(source)=='534a0f5bcd676af37cda6af9b7085da012dc68edc90ac99c3fdffc06b971e341'
binary=Path('/Users/n/numivivo-duration-interpolation-20260912/temporal')
a=np.load(source);training=[]
for d in ['D34','D38','D39']:
 prefix='query:GSE226572:'+d;c=a[prefix+'|noChange']
 for key in a.files:
  if key.startswith(prefix+'|truth|IFNB:'):training.append({'donor':d,'time':float(key.rsplit(':',1)[1][:-1]),'response':(a[key]-c).tolist()})
k=anndata.read_h5ad(r/'training.h5ad');ids=k.var_names.tolist();indices=[ids.index(str(x)) for x in a['featureIDs']];samples={x['id']:x for x in json.loads((r/'training.json').read_text())['mapping']['samples']}
# Normalize with the full source feature denominator before projection.
x=k.X.toarray() if hasattr(k.X,'toarray') else np.asarray(k.X);v=np.log1p(x/x.sum(1)[:,None]*1e6)[:,indices];donors=sorted({x['donorID'] for x in samples.values()});assert len(donors)==8
freeze=json.loads((r/'freeze.json').read_text())
for i,record in enumerate(freeze['records']):
 assert h(r/'inputs'/f'{i}.json')==record['inputSHA256']
 assert h(r/'outputs'/f'{i}.json')==record['outputSHA256']
curve={t:np.mean([z['response'] for z in training if z['time']==t],axis=0) for t in sorted({z['time'] for z in training})};response=curve[4]*.5+curve[8]*.5
baselineResponse=np.mean([np.mean([z['response'] for z in training if z['donor']==d],axis=0) for d in sorted({z['donor'] for z in training})],axis=0)
rows=[]
for i,d in enumerate(donors):
 ci=[j for j,s in enumerate(k.obs['sample']) if samples[s]['donorID']==d and samples[s]['condition']=='control'][0];ti=[j for j,s in enumerate(k.obs['sample']) if samples[s]['donorID']==d and samples[s]['condition']=='IFNB'];assert len(ti)==1;b=json.loads((r/'outputs'/f'{i}.json').read_text());err=float(np.max(abs(np.maximum(0,v[ci]+response)-b['predicted'])));assert err<1e-12
 baselineError=float(np.max(abs(np.maximum(0,v[ci]+baselineResponse)-b['baseline'])));assert baselineError<1e-12
 scores={key:math.sqrt(math.fsum((float(a)-float(c))**2 for a,c in zip(pred,v[ti[0]]))/len(indices)) for key,pred in [('candidate',b['predicted']),('trainingMean',b['baseline']),('noChange',v[ci])]};gains={key:1-scores['candidate']/scores[key] for key in ['trainingMean','noChange']};rows.append({'donor':d,'rmse':scores,'gains':gains,'pass':all(x>=.05 for x in gains.values()),'maxReferenceDifference':err,'maxBaselineReferenceDifference':baselineError})
(r/'results.json').write_text(json.dumps({'rows':rows,'passingDonors':sum(x['pass'] for x in rows),'allDonorsPass':all(x['pass'] for x in rows)},indent=2));print([(x['donor'],x['gains'],x['pass']) for x in rows])
