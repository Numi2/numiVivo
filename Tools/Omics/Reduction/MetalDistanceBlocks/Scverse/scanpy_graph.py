from pathlib import Path
import os,time,sys,json
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',NUMBA_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import numpy as np,scanpy as sc,anndata as ad
from scipy import sparse
from threadpoolctl import threadpool_limits
source=Path(sys.argv[1]);out=Path(sys.argv[2]);sc.settings.n_jobs=1
with threadpool_limits(limits=1):
 start=time.perf_counter();a=np.fromfile(source,dtype=[('row','<u4'),('col','<u4'),('value','<f8')]);n=24673;d=20
 assert len(a)==n*d and np.array_equal(a['row'],np.repeat(np.arange(n),d)) and np.array_equal(a['col'],np.tile(np.arange(d),n))
 x=a['value'].reshape(n,d).copy();data=ad.AnnData(X=sparse.csr_matrix((n,0)),obsm={'X_pca':x})
 sc.pp.neighbors(data,n_neighbors=20,n_pcs=20,use_rep='X_pca',knn=True,method='umap',metric='euclidean',transformer='sklearn',random_state=7)
 elapsed=time.perf_counter()-start
 out.mkdir();sparse.save_npz(out/'distances.npz',data.obsp['distances']);sparse.save_npz(out/'connectivities.npz',data.obsp['connectivities'])
 (out/'timing.json').write_text(json.dumps(dict(backend='scanpy',graphStageSeconds=elapsed,cells=n,sourceScoreReadsIncluded=True,params=data.uns['neighbors']['params']),indent=2))
 print(elapsed)
