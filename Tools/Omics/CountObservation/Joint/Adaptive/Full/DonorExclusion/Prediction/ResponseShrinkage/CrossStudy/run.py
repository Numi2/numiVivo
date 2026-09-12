from pathlib import Path
import json,hashlib,subprocess,time
import anndata as ad,numpy as np
s=Path(__file__).parent;source=Path('/Users/home/numivivo-cross-study-ifnb-20260911/inputs');freeze=json.loads((source/'input-freeze.json').read_text());panel=json.loads((source/'panel.json').read_text());folds=json.loads((source/'folds.json').read_text());out=s/'native';out.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bound(name):
 p=source/name;assert sha(p)==freeze['files'][name],name;return p
for n in ['panel.json','folds.json']:bound(n)
record={'inputFreezeSHA256':sha(source/'input-freeze.json'),'executableSHA256':sha(s/'response-shrinkage'),'scoringStarted':False,'folds':{}}
for fold in folds:
 tag=fold['id'];a=ad.read_h5ad(bound(tag+'/training.h5ad'));q=ad.read_h5ad(bound(tag+'/query.h5ad'));plan=json.loads(bound(tag+'/training.json').read_text());qp=json.loads(bound(tag+'/query.json').read_text());samples={v['id']:v for v in plan['mapping']['samples']};assert len(samples)==a.n_obs
 lookup={(samples[str(v)]['donorID'],samples[str(v)]['condition']):i for i,v in enumerate(a.obs['sample'])};donors=sorted(set(d for d,c in lookup));assert len(lookup)==2*len(donors)
 assert fold['heldOutDonor'] not in donors and qp['mapping']['samples'][0]['donorID']==fold['heldOutDonor'] and qp['mapping']['samples'][0]['condition']=='control' and q.n_obs==1
 def logs(a):
  values=a.X.toarray() if hasattr(a.X,'toarray') else np.asarray(a.X);assert np.all(values>=0);totals=values.sum(axis=1,dtype=np.uint64);assert np.all(totals>0);normalized=np.log1p(values.astype('f8')/totals[:,None]*1e6);indices=a.var_names.get_indexer(panel);assert np.all(indices>=0) and a.var_names.is_unique;return normalized[:,indices]
 train=logs(a);query=logs(q)[0];x=train[[lookup[(d,'control')] for d in donors]];y=train[[lookup[(d,'IFNB')] for d in donors]]
 obj={'featureIDs':panel,'trainingDonorIDs':donors,'controls':x.tolist(),'treated':y.tolist(),'queryDonorID':fold['heldOutDonor'],'queryControl':query.tolist()};inp=out/(tag+'-input.json');dest=out/(tag+'-output.json');inp.write_text(json.dumps(obj,separators=(',',':'))+'\n');start=time.time();subprocess.run([str(s/'response-shrinkage'),str(inp),str(dest)],check=True)
 record['folds'][tag]={'inputSHA256':sha(inp),'outputSHA256':sha(dest),'seconds':time.time()-start};print(tag,flush=True)
record['status']='complete';(s/'prediction-freeze.json').write_text(json.dumps(record,indent=2)+'\n')
