from pathlib import Path
import numpy as np, h5py,json,hashlib
from scipy import sparse
s=Path('/Users/n/numivivo-streamed-logcpm-20260912')
source=Path('/Users/n/numivivo-paired-workflow-20260912/bundle/original.h5')
assert hashlib.sha256(source.read_bytes()).hexdigest()=='5fbff5a4d85e0df345f6502e966ec787a8a4c429fd6b88a8772c43fd915cf3ff'
with h5py.File(source) as f:
 g=f['matrix']; x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:])).T.tocsr()
 mask=g['features/feature_type'].asstr()[:]=='Gene Expression'
 ids=g['features/id'].asstr()[:][mask].tolist();x=x[:,mask];x.sort_indices()
assert x.shape==(2711,36601) and x.nnz==5218473 and x.sum()==11786194
totals=np.asarray(x.sum(axis=1)).ravel();groups=np.arange(x.shape[0])%3
(s/'plan.json').write_text(json.dumps(dict(featureIDs=ids,groupIDs=['partition0','partition1','partition2'],rowGroups=groups.tolist(),rowTotals=totals.tolist())))
with (s/'counts.bin').open('xb') as f:
 for i in range(x.shape[0]):
  lo,hi=x.indptr[i:i+2];a=np.empty(hi-lo,dtype=[('row','<u4'),('feature','<u4'),('count','<u8')]);a['row']=i;a['feature']=x.indices[lo:hi];a['count']=x.data[lo:hi];f.write(a.tobytes())
y=x.astype(float);y.data=np.log1p(y.data/np.repeat(totals,np.diff(y.indptr))*1e6)
ref=np.stack([np.asarray(y[groups==j].mean(axis=0)).ravel() for j in range(3)])
np.save(s/'reference.npy',ref)
print('PASS prepared full RNA source; partitions are numerical checks only')
