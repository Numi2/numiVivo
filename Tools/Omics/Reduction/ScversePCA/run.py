from pathlib import Path
import json,time,hashlib,inspect,importlib.metadata
import numpy as np, scipy, scanpy as sc, anndata
r=Path('/Users/n/numivivo-pca-scverse-20260912')
protocol={'scope':'Full cohort sparse normalized matrix and same native-selected 2000 features; numerical PCA comparison, not matched full-pipeline timing','cohorts':['kang','hagai'],'components':20,'solver':'arpack','seed':7,'dtype':'float64','minimumSubspaceSingularValue':0.99999,'maximumRelativeVarianceDifference':1e-5}
(r/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
(r/'scanpy-pca-source.py').write_text(inspect.getsource(sc.pp.pca))
(r/'versions.json').write_text(json.dumps({x:importlib.metadata.version(x) for x in ['scanpy','anndata','numpy','scipy','scikit-learn']},indent=2)+'\n')
def h(p):
 s=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):s.update(b)
 return s.hexdigest()
def records(path,rows,cols):
 a=np.fromfile(path,dtype=[('row','<u4'),('col','<u4'),('value','<f8')]);assert len(a)==rows*cols
 assert np.array_equal(a['row'],np.repeat(np.arange(rows),cols)) and np.array_equal(a['col'],np.tile(np.arange(cols),rows))
 return a['value'].reshape(rows,cols)
results=[]
for cohort in protocol['cohorts']:
 p=Path('/Users/n/numivivo-pca-borrowed-20260912')/('pca-1' if cohort=='kang' else 'hagai-pca-1')
 model=json.loads((p/'model.json').read_text());metadata=json.loads((p/'metadata.json').read_text())
 start=time.monotonic();a=anndata.read_h5ad(p/'original.h5ad');assert scipy.sparse.issparse(a.X)
 mapping=json.loads((p/'plan.json').read_text())['mapping']
 barcodes=a.obs[mapping['barcodeColumn']].astype(str).tolist() if 'barcodeColumn' in mapping else a.obs_names.tolist()
 assert barcodes==[x['barcode'] for x in metadata['cells']]
 assert a.obs[mapping['sampleColumn']].astype(str).tolist()==[x['sampleID'] for x in metadata['cells']]
 a.X=a.X.astype(np.float64);sc.pp.normalize_total(a,target_sum=model['normalizationTarget']);sc.pp.log1p(a)
 selected=model['selectedFeatureIndices'];x=a.X[:,selected].tocsr();preparation=time.monotonic()-start
 means=np.asarray(x.mean(axis=0)).ravel();center_error=float(np.max(np.abs(means-np.asarray(model['projectionCenters']))))
 start=time.monotonic();scores,loadings,ratio,variance=sc.pp.pca(x,n_comps=20,zero_center=True,svd_solver='arpack',random_state=7,return_info=True,dtype='float64');elapsed=time.monotonic()-start
 native=records(p/'loadings.bin',len(selected),20);native_scores=records(p/'scores.bin',len(a.obs),20)
 ref=loadings.T;singular=np.linalg.svd(native.T@ref,compute_uv=False)
 signs=np.sign(np.sum(native*ref,axis=0));aligned=ref*signs
 vdiff=float(np.max(np.abs(variance/np.array(model['explainedVariance'])-1)))
 item={'cohort':cohort,'cells':len(a.obs),'selectedFeatures':len(selected),'sourceSHA256':h(p/'original.h5ad'),'modelSHA256':h(p/'model.json'),'loadingsSHA256':h(p/'loadings.bin'),'scoresSHA256':h(p/'scores.bin'),'preparationSeconds':preparation,'scanpyPCASeconds':elapsed,'minimumSubspaceSingularValue':float(singular.min()),'maximumRelativeVarianceDifference':vdiff,'maximumCenterDifference':center_error,'maximumSignAlignedLoadingDifference':float(np.max(np.abs(native-aligned))),'maximumSignAlignedScoreDifference':float(np.max(np.abs(native_scores-scores*signs)))}
 item['gatePass']=item['minimumSubspaceSingularValue']>=protocol['minimumSubspaceSingularValue'] and vdiff<=protocol['maximumRelativeVarianceDifference'];results.append(item)
 np.savez_compressed(r/(cohort+'-reference.npz'),scores=scores,loadings=loadings,variance=variance,ratio=ratio)
 (r/'results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(item),flush=True)
