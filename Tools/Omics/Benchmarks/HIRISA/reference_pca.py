#!/usr/bin/env python3
"""Complete-source HIRISA centered PCA with a file-backed SciPy sparse operator.

This independent reference never loads the raw cells-by-genes matrix or native
PCA outputs. Its selected sparse values retain every row and normalize using all
original genes. Numerical results are frozen before native comparison.
"""
import argparse, hashlib, importlib.metadata, json, os, platform, time
from pathlib import Path
import h5py
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator, eigsh

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''): h.update(block)
 return h.hexdigest()

def write(path,value):
 with path.open('x') as f: json.dump(value,f,indent=2,sort_keys=True,allow_nan=False); f.write('\n')

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument('--root',type=Path,required=True)
 a=p.parse_args(); root=a.root; measured=root/'pca'; out=measured/'independent-reference'
 out.mkdir(exist_ok=False); start=time.monotonic(); source=root/'hirisa.h5ad'
 measurement=json.loads((measured/'measurement.json').read_text())
 plan=json.loads((measured/'full-plan.json').read_text()); settings=plan['reduction']['pca']
 assert sha(source)==measurement['sourceSHA256']
 assert sha(measured/'feature-moments.npz')==measurement['featureMomentsSHA256']
 moments=np.load(measured/'feature-moments.npz',allow_pickle=False)
 selected=moments['selectedFeatureIndices']; n=measurement['cells']; width=len(selected); nnz=measurement['selectedEntries']
 spec=dict(schemaVersion=1,sourceSHA256=measurement['sourceSHA256'],measurementSHA256=sha(measured/'measurement.json'),
  planSHA256=sha(measured/'full-plan.json'),scriptSHA256=sha(Path(__file__)),cells=n,features=measurement['features'],
  selectedFeatures=width,selectedEntries=nnz,solver='scipy.sparse.linalg.eigsh centered covariance LinearOperator',
  components=settings['components'],seed=7,tolerance=1e-10,maximumIterations=2000,
  chunkRows=1024,host=platform.node(),threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
  versions={k:importlib.metadata.version(k) for k in ('numpy','scipy','h5py')},
  scope='Complete-source numerical reference; no native outputs read, no dense cells by genes matrix, no biological or controlled cross-host performance claim.')
 write(out/'specification.json',spec)
 data=np.lib.format.open_memmap(out/'values.npy',mode='w+',dtype='<f8',shape=(nnz,))
 indices=np.lib.format.open_memmap(out/'indices.npy',mode='w+',dtype='<i4',shape=(nnz,))
 indptr=np.lib.format.open_memmap(out/'indptr.npy',mode='w+',dtype='<i4',shape=(n+1,)); indptr[0]=0
 assert nnz<np.iinfo(np.int32).max
 lookup=np.full(measurement['features'],-1,dtype=np.int32); lookup[selected]=np.arange(width,dtype=np.int32)
 offset=0; raw_entries=0; umis=0; sums=np.zeros(width); squares=np.zeros(width)
 with h5py.File(source,'r') as f:
  x=f['X']; assert tuple(x.attrs['shape'])==(n,measurement['features'])
  for first in range(0,n,spec['chunkRows']):
   last=min(n,first+spec['chunkRows']); ptr=x['indptr'][first:last+1]; left,right=int(ptr[0]),int(ptr[-1]); ptr=ptr-ptr[0]
   raw=x['data'][left:right]; cols=x['indices'][left:right]; rows=np.repeat(np.arange(last-first),np.diff(ptr))
   totals=np.bincount(rows,weights=raw,minlength=last-first)
   np.testing.assert_array_equal(totals,f['obs/n_umis'][first:last])
   local=lookup[cols]; keep=local>=0; retained_rows=rows[keep]; retained_cols=local[keep]
   values=np.log1p(raw[keep].astype(np.float64)/totals[retained_rows]*plan['reduction']['normalizationTarget'])
   end=offset+len(values); assert end<=nnz
   data[offset:end]=values; indices[offset:end]=retained_cols
   indptr[first+1:last+1]=offset+np.cumsum(np.bincount(retained_rows,minlength=last-first))
   sums+=np.bincount(retained_cols,weights=values,minlength=width)
   squares+=np.bincount(retained_cols,weights=values*values,minlength=width)
   offset=end; raw_entries+=len(raw); umis+=int(totals.sum())
   if last%131072==0 or last==n: print('materializedCells',last,'selectedEntries',offset,flush=True)
 assert offset==nnz and raw_entries==measurement['canonicalNonzeros'] and umis==measurement['UMIs']
 for array in (data,indices,indptr): array.flush()
 centers=sums/n; variances=(squares-sums*sums/n)/(n-1)
 np.testing.assert_allclose(centers,moments['meanLogNormalized'][selected],rtol=1e-11,atol=1e-13)
 np.testing.assert_allclose(variances,moments['varianceLogNormalized'][selected],rtol=1e-10,atol=1e-13)
 matrix=csr_matrix((data,indices,indptr),shape=(n,width),copy=False)
 assert np.shares_memory(matrix.data,data) and np.shares_memory(matrix.indices,indices) and np.shares_memory(matrix.indptr,indptr)
 assert matrix.has_sorted_indices and matrix.has_canonical_format
 materialized=time.monotonic(); calls=0
 def covariance(v):
  nonlocal calls
  projected=matrix@v-float(centers@v)
  result=(matrix.T@projected-centers*projected.sum())/(n-1)
  calls+=1
  if calls%20==0: print('covarianceApplications',calls,'seconds',time.monotonic()-materialized,flush=True)
  return result
 operator=LinearOperator((width,width),matvec=covariance,dtype=np.float64)
 eigenvalues,loadings=eigsh(operator,k=spec['components'],which='LA',tol=spec['tolerance'],maxiter=spec['maximumIterations'],v0=np.random.default_rng(spec['seed']).normal(size=width))
 order=np.argsort(eigenvalues)[::-1]; eigenvalues=eigenvalues[order]; loadings=loadings[:,order]
 for j in range(len(eigenvalues)):
  pivot=np.argmax(np.abs(loadings[:,j])); loadings[:,j]*=1 if loadings[pivot,j]>=0 else -1
 residuals=np.array([np.linalg.norm(covariance(loadings[:,j])-eigenvalues[j]*loadings[:,j])/max(abs(eigenvalues[j]),1e-300) for j in range(len(eigenvalues))])
 orthogonality=float(np.max(np.abs(loadings.T@loadings-np.eye(len(eigenvalues)))))
 assert np.isfinite(eigenvalues).all() and (eigenvalues>0).all() and residuals.max()<1e-8 and orthogonality<1e-10
 np.savez_compressed(out/'fit.npz',selectedFeatureIndices=selected,projectionCenters=centers,varianceLogNormalized=variances,
  explainedVariance=eigenvalues,explainedVarianceRatio=eigenvalues/variances.sum(),loadings=loadings,relativeResiduals=residuals)
 result=dict(schemaVersion=1,status='passed-reference-fit',specificationSHA256=sha(out/'specification.json'),fitSHA256=sha(out/'fit.npz'),
  materializationSeconds=materialized-start,totalSeconds=time.monotonic()-start,covarianceApplications=calls,
  maximumRelativeResidual=float(residuals.max()),maximumLoadingOrthogonalityError=orthogonality,
  scratch={p.name:dict(bytes=p.stat().st_size,SHA256=sha(p)) for p in (out/'values.npy',out/'indices.npy',out/'indptr.npy')},
  nativeComparisonPerformed=False)
 write(out/'result.json',result); print(json.dumps(result),flush=True)

if __name__=='__main__': main()
