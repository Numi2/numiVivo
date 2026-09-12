from pathlib import Path
import json,hashlib,subprocess,math
import numpy as np,anndata
r=Path(__file__).parent;h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();source=Path('/Users/n/numivivo-duration-source-20260912.npz');assert h(source)=='534a0f5bcd676af37cda6af9b7085da012dc68edc90ac99c3fdffc06b971e341'
binary=Path('/Users/n/numivivo-duration-interpolation-20260912/temporal')
(r/'protocol.json').write_text(json.dumps({'scope':'Reused development cross-study whole-PBMC transfer GSE226572 to Kang','time':6,'gate':'at least 5 percent lower RMSE than both mean response and no change in each donor','sourceSHA256':h(source),'kangH5ADSHA256':h(r/'training.h5ad'),'kangMappingSHA256':h(r/'training.json'),'binarySHA256':h(binary),'model':'unchanged temporal interpolation; no target outcome selection or refit'},indent=2))
a=np.load(source);training=[]
for d in ['D34','D38','D39']:
 prefix='query:GSE226572:'+d;c=a[prefix+'|noChange']
 for key in a.files:
  if key.startswith(prefix+'|truth|IFNB:'):training.append({'donor':d,'time':float(key.rsplit(':',1)[1][:-1]),'response':(a[key]-c).tolist()})
k=anndata.read_h5ad(r/'training.h5ad');ids=k.var_names.tolist();indices=[ids.index(str(x)) for x in a['featureIDs']];samples={x['id']:x for x in json.loads((r/'training.json').read_text())['mapping']['samples']}
# Normalize with the full source feature denominator before projection.
x=k.X.toarray() if hasattr(k.X,'toarray') else np.asarray(k.X);v=np.log1p(x/x.sum(1)[:,None]*1e6)[:,indices];donors=sorted({x['donorID'] for x in samples.values()});assert len(donors)==8
(r/'inputs').mkdir();(r/'outputs').mkdir();records=[]
for i,d in enumerate(donors):
 ix=[j for j,s in enumerate(k.obs['sample']) if samples[s]['donorID']==d and samples[s]['condition']=='control'];assert len(ix)==1
 inp=r/'inputs'/f'{i}.json';out=r/'outputs'/f'{i}.json';inp.write_text(json.dumps({'donor':d,'time':6,'control':v[ix[0]].tolist(),'training':training}));subprocess.run([str(binary),str(inp),str(out)],check=True);records.append({'donor':d,'inputSHA256':h(inp),'outputSHA256':h(out)})
(r/'freeze.json').write_text(json.dumps({'scored':False,'records':records},indent=2))
curve={t:np.mean([z['response'] for z in training if z['time']==t],axis=0) for t in sorted({z['time'] for z in training})};response=curve[4]*.75+curve[8]*.25
rows=[]
for i,d in enumerate(donors):
 ci=[j for j,s in enumerate(k.obs['sample']) if samples[s]['donorID']==d and samples[s]['condition']=='control'][0];ti=[j for j,s in enumerate(k.obs['sample']) if samples[s]['donorID']==d and samples[s]['condition']=='IFNB'];assert len(ti)==1;b=json.loads((r/'outputs'/f'{i}.json').read_text());err=float(np.max(abs(np.maximum(0,v[ci]+response)-b['predicted'])));assert err<1e-12
 scores={key:math.sqrt(math.fsum((float(a)-float(c))**2 for a,c in zip(pred,v[ti[0]]))/len(indices)) for key,pred in [('candidate',b['predicted']),('trainingMean',b['baseline']),('noChange',v[ci])]};gains={key:1-scores['candidate']/scores[key] for key in ['trainingMean','noChange']};rows.append({'donor':d,'rmse':scores,'gains':gains,'pass':all(x>=.05 for x in gains.values()),'maxReferenceDifference':err})
(r/'results.json').write_text(json.dumps({'rows':rows,'passingDonors':sum(x['pass'] for x in rows),'allDonorsPass':all(x['pass'] for x in rows)},indent=2));print([(x['donor'],x['gains'],x['pass']) for x in rows])
