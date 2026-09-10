#!/usr/bin/env python3
"""Compare every native HIRISA PCA score against the frozen sparse reference."""
import argparse, hashlib, itertools, json, time
from pathlib import Path
import h5py
import ijson
import numpy as np
from anndata.io import read_elem
from scipy.sparse import csr_matrix

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
 a=p.parse_args();assert not a.out.exists();start=time.monotonic();base=a.root/'pca';ref=base/'independent-reference';bundle=a.bundle
 measurement=json.loads((base/'measurement.json').read_text());prior=json.loads((base/'native-json-measurement.json').read_text())
 fit_receipt=json.loads((ref/'result.json').read_text());assert fit_receipt['status']=='passed-reference-fit'
 assert sha(ref/'fit.npz')==fit_receipt['fitSHA256'];assert sha(ref/'specification.json')==fit_receipt['specificationSHA256']
 receipt=json.loads((bundle/'receipt.json').read_text());assert receipt['matrixFormat']=='complete-row-major-u32-row-u32-component-f64-le/v1'
 assert bytes(receipt['source']['bytes']).hex()==measurement['sourceSHA256']
 artifact_hashes={}
 for key in ('plan','metadata','quality','model','scores','loadings'):
  path=bundle/(key+('.bin' if key in ('scores','loadings') else '.json'));artifact_hashes[key]=sha(path)
  assert artifact_hashes[key]==bytes(receipt[key]['bytes']).hex(),key
 for key in ('metadata','quality'):
  assert artifact_hashes[key]==prior[key]['SHA256'] and (bundle/(key+'.json')).stat().st_size==prior[key]['bytes']
 model=json.loads((bundle/'model.json').read_text());plan=json.loads((bundle/'plan.json').read_text())
 assert plan==json.loads((base/'full-plan.json').read_text())
 fit=np.load(ref/'fit.npz',allow_pickle=False);moments=np.load(base/'feature-moments.npz',allow_pickle=False)
 n=measurement['cells'];selected=fit['selectedFeatureIndices'];k=model['options']['components'];width=len(selected)
 assert model['cells']==n and model['canonicalNonzeros']==measurement['canonicalNonzeros']
 np.testing.assert_array_equal(model['selectedFeatureIndices'],selected)
 features=model['features'];assert len(features)==measurement['features']
 assert [f['featureID'] for f in features]==moments['featureIDs'].tolist()
 assert [f['featureIndex'] for f in features]==list(range(len(features)))
 np.testing.assert_array_equal([f['selected'] for f in features],np.isin(np.arange(len(features)),selected))
 np.testing.assert_array_equal([f['meanBin'] for f in features],moments['meanBins'])
 errors={}
 for key in ('meanNormalized','varianceNormalized','logMeanForBinning','logDispersion','normalizedDispersion'):
  observed=np.array([np.nan if f.get(key) is None else f[key] for f in features]);expected=moments[key]
  np.testing.assert_array_equal(np.isnan(observed),np.isnan(expected))
  tolerance=2e-6 if key=='normalizedDispersion' else 1e-9
  np.testing.assert_allclose(observed,expected,rtol=tolerance,atol=tolerance,equal_nan=True)
  errors[key]=float(np.nanmax(np.abs(observed-expected)))
 centers=np.asarray(model['projectionCenters']);np.testing.assert_allclose(centers,fit['projectionCenters'],rtol=1e-12,atol=1e-12)
 variance=np.asarray(model['explainedVariance']);np.testing.assert_allclose(variance,fit['explainedVariance'],rtol=1e-7,atol=1e-9)
 np.testing.assert_allclose(model['explainedVarianceRatio'],fit['explainedVarianceRatio'],rtol=1e-7,atol=1e-10)
 dtype=np.dtype([('row','<u4'),('component','<u4'),('value','<f8')])
 def records(name,rows):
  file=bundle/name;assert file.stat().st_size==rows*k*16
  return np.memmap(file,mode='r',dtype=dtype,shape=(rows,k))
 def values(records,first,last):
  block=records[first:last];np.testing.assert_array_equal(block['row'],np.broadcast_to(np.arange(first,last)[:,None],(last-first,k)))
  np.testing.assert_array_equal(block['component'],np.broadcast_to(np.arange(k),(last-first,k)))
  v=block['value'].copy();assert np.isfinite(v).all();return v
 loading=values(records('loadings.bin',width),0,width);reference_loading=fit['loadings']
 u,cosines,vt=np.linalg.svd(loading.T@reference_loading);rotation=u@vt
 assert cosines.min()>1-1e-8
 orthogonality=float(np.abs(loading.T@loading-np.eye(k)).max());assert orthogonality<1e-8
 arrays={name:np.load(ref/name,mmap_mode='r',allow_pickle=False) for name in ('values.npy','indices.npy','indptr.npy')}
 for name,array in arrays.items(): assert sha(ref/name)==fit_receipt['scratch'][name]['SHA256']
 matrix=csr_matrix((arrays['values.npy'],arrays['indices.npy'],arrays['indptr.npy']),shape=(n,width),copy=False)
 assert np.shares_memory(matrix.data,arrays['values.npy']) and np.shares_memory(matrix.indices,arrays['indices.npy'])
 def covariance(v):
  projected=matrix@v-float(fit['projectionCenters']@v)
  return (matrix.T@projected-fit['projectionCenters']*projected.sum())/(n-1)
 residuals=np.array([np.linalg.norm(covariance(loading[:,j])-variance[j]*loading[:,j])/variance[j] for j in range(k)])
 assert residuals.max()<=model['options']['relativeResidualTolerance']
 assert max(model['relativeResiduals'])<=model['options']['relativeResidualTolerance']
 score_records=records('scores.bin',n);sums=np.zeros(k);cross=np.zeros((k,k));error2=0.0;norm2=0.0;max_projection_error=0.0
 source=a.root/'hirisa.h5ad';assert sha(source)==measurement['sourceSHA256']
 with h5py.File(source,'r') as h, (bundle/'metadata.json').open('rb') as mf, (bundle/'quality.json').open('rb') as qf:
  obs=h['obs'];cells=ijson.items(mf,'cells.item',use_float=True);quality=ijson.items(qf,'item',use_float=True)
  for first in range(0,n,4096):
   last=min(n,first+4096);count=last-first;native=values(score_records,first,last)
   block=matrix[first:last];expected=block@reference_loading-fit['projectionCenters']@reference_loading
   delta=native@rotation-expected;error2+=float(np.sum(delta*delta));norm2+=float(np.sum(expected*expected))
   direct=block@loading-centers@loading;max_projection_error=max(max_projection_error,float(np.abs(native-direct).max()))
   np.testing.assert_allclose(native,direct,rtol=1e-9,atol=1e-9)
   sums+=native.sum(axis=0);cross+=native.T@native
   ids=obs['_index'].asstr()[first:last];samples=obs['geo_accession'].asstr()[first:last]
   qc={key:obs[key][first:last] for key in ('n_umis','n_genes','n_mito_umis')}
   meta=list(itertools.islice(cells,count));qual=list(itertools.islice(quality,count));assert len(meta)==len(qual)==count
   for j,(c,q) in enumerate(zip(meta,qual,strict=True)):
    assert c['barcode']==q['barcode']==ids[j] and c['sampleID']==q['sampleID']==samples[j]
    assert q['totalCounts']==int(qc['n_umis'][j]) and q['detectedFeatures']==int(qc['n_genes'][j]) and q['mitochondrialCounts']==int(qc['n_mito_umis'][j])
    assert q['mitochondrialFeatureCount']==11
   if last%262144==0 or last==n:print('comparedRows',last,flush=True)
  sentinel=object();assert next(cells,sentinel) is sentinel and next(quality,sentinel) is sentinel
 score_error=float(np.sqrt(error2/norm2));assert score_error<1e-5
 np.testing.assert_allclose(sums/n,0,atol=1e-10)
 np.testing.assert_allclose(cross/(n-1),np.diag(variance),rtol=1e-7,atol=1e-8)
 storage=model['storage'];assert storage['sourcePasses']==3
 assert storage['selectedEntries']==measurement['selectedEntries'] and storage['cacheBytes']==measurement['currentCOOCacheBytes']
 assert storage['entryVisits']==measurement['plannedEntryVisits']
 result=dict(schemaVersion=1,status='passed',cells=n,features=len(features),selectedFeatures=width,components=k,
  sourceSHA256=measurement['sourceSHA256'],artifacts=artifact_hashes,referenceFitSHA256=fit_receipt['fitSHA256'],
  comparisonSpecificationSHA256=sha(Path(__file__).with_name('PCA_COMPARISON.md')),comparisonScriptSHA256=sha(Path(__file__)),
  featureMaximumAbsoluteErrors=errors,minimumSubspaceCosine=float(cosines.min()),alignedRelativeScoreError=score_error,
  maximumDirectProjectionAbsoluteError=max_projection_error,maximumNativeResidual=max(model['relativeResiduals']),
  maximumIndependentlyMeasuredNativeResidual=float(residuals.max()),maximumLoadingOrthogonalityError=orthogonality,
  allOriginalIdentitiesAndQualityChecked=True,allScoreCoordinatesAndValuesChecked=True,exactPriorNativeMetadataAndQuality=True,
  storage=storage,seconds=time.monotonic()-start,
  scope='Complete-source CPU numerical qualification only. No graph, integration, biological validity, GPU or controlled cross-host performance claim.')
 with a.out.open('x') as f:json.dump(result,f,sort_keys=True,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps(result),flush=True)

if __name__=='__main__':main()
