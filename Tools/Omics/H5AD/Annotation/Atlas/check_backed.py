from pathlib import Path
import argparse,json,time,warnings,gc
from importlib.metadata import version
import anndata as ad,h5py,numpy as np,pandas as pd
from scipy import sparse
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--annotated',type=Path,required=True);p.add_argument('--report',type=Path,required=True);args=p.parse_args();start=time.time();checks=[]
with warnings.catch_warnings():warnings.simplefilter('ignore');atlas=ad.read_h5ad(args.annotated,backed='r')
try:
 assert atlas.isbacked and getattr(atlas.X,'format',None)=='csr' and not sparse.issparse(atlas.X)
 n,m=atlas.shape;print('opened',n,m,type(atlas.X).__name__,flush=True)
 np.testing.assert_array_equal(atlas.obs['numivivo_source_row'].to_numpy(),np.arange(n,dtype='i8'))
 with h5py.File(args.source,'r') as f:
  for frame_name,frame in [('obs',atlas.obs),('var',atlas.var)]:
   g=f[frame_name];key=g.attrs['_index'];index=pd.Index(ad.io.read_elem(g[key]),name=None if key=='_index' else key);pd.testing.assert_index_equal(frame.index,index);del index
   for column in g.attrs['column-order']:
    expected=pd.Series(ad.io.read_elem(g[column]),index=frame.index,name=column);pd.testing.assert_series_equal(frame[column],expected);checks.append(frame_name+'/'+column);del expected;gc.collect();print('column',frame_name,column,flush=True)
  selected=sorted({0,1,7,n//2,n-2,n-1});matrix=f['X'];count=0
  for row in selected:
   lo,hi=map(int,matrix['indptr'][row:row+2]);expected=sparse.csr_matrix((matrix['data'][lo:hi],matrix['indices'][lo:hi],np.array([0,hi-lo])),shape=(1,m));actual=atlas.X[row:row+1,:];assert actual.shape==expected.shape
   np.testing.assert_array_equal(actual.indptr,expected.indptr);np.testing.assert_array_equal(actual.indices,expected.indices);np.testing.assert_array_equal(actual.data,expected.data);count+=hi-lo
 report={'status':'passed-complete-backed-open-and-annotations','source':str(args.source),'annotated':str(args.annotated),'cells':n,'features':m,'matrixBackend':type(atlas.X).__module__+'.'+type(atlas.X).__name__,'allOriginalObservationAndFeatureColumnsChecked':checks,'newRowColumnExact':True,'sampledCountRows':selected,'sampledCountEntries':count,'anndata':version('anndata'),'pandas':version('pandas'),'h5py':version('h5py'),'seconds':time.time()-start,'scope':'Complete backed AnnData open and full obs/var metadata equivalence; count slicing is sampled here. The separate HDF5 verifier checked every original stored element. AnnData metadata remain resident; not a biological prediction or end-to-end speedup test.'};args.report.write_text(json.dumps(report,indent=2)+'\n');print('PASS',flush=True)
finally:atlas.file.close()
