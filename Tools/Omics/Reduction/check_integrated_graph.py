#!/usr/bin/env python3
"""All-query exact-neighbor recall and independent integrated graph invariants."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from scipy import sparse
from scipy.sparse.csgraph import connected_components
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
root=a.root/'graph';g=json.loads((root/'graph.json').read_text());assert g['options']['representation']=='integrated';n=g['cells'];k=g['options']['neighbors'];d=g['dimensions'];ids=np.arange(n)
def records(path):return np.fromfile(path,dtype=[('row','<u4'),('column','<u4'),('value','<f8')])
r=records(root/'input/scores.bin').reshape(n,d);assert np.equal(r['row'],ids[:,None]).all() and np.equal(r['column'],np.arange(d)).all();x=r['value'].copy();neighbors=records(root/'neighbors.bin').reshape(n,k)
assert np.equal(neighbors['row'],ids[:,None]).all() and np.equal(neighbors['column'][:,0],ids).all() and np.equal(neighbors['value'][:,0],0).all()
recalls=[];maximum_error=0
for start in range(0,n,64):
 distance=cdist(x[start:start+64],x)
 for off,row in enumerate(distance):
  i=start+off;row[i]=np.inf;threshold=np.partition(row,k-2)[k-2];candidates=ids[row<=threshold];exact=candidates[np.lexsort((candidates,row[candidates]))[:k-1]];actual=neighbors['column'][i,1:]
  assert len(set(actual))==k-1 and (actual!=i).all() and (actual<n).all()
  recalls.append(len(np.intersect1d(actual,exact))/(k-1));delta=float(np.max(np.abs(neighbors['value'][i,1:]-row[actual])));maximum_error=max(maximum_error,delta)
  np.testing.assert_allclose(neighbors['value'][i,1:],row[actual],rtol=1e-12,atol=1e-12)
np.save(a.out/'per-cell-recall.npy',recalls)
edges=records(root/'edges.bin');adj=sparse.csr_matrix((edges['value'],(edges['row'],edges['column'])),shape=(n,n));assert np.isfinite(adj.data).all() and (adj.data>0).all() and (adj.diagonal()==0).all() and (adj!=adj.T).nnz==0
components=connected_components(adj,directed=False)[0];assert components==g['connectedComponents']
result=dict(status='passed' if np.mean(recalls)>=0.95 and np.percentile(recalls,5)>=0.85 else 'recall-gap',cells=n,queryCells=n,neighborsExcludingSelf=k-1,meanRecall=float(np.mean(recalls)),p5Recall=float(np.percentile(recalls,5)),maximumDistanceError=maximum_error,connectedComponents=components,allPositiveFiniteSymmetricEdges=True,graphSHA256=hashlib.sha256((root/'graph.json').read_bytes()).hexdigest(),scoreSHA256=hashlib.sha256((root/'input/scores.bin').read_bytes()).hexdigest(),qualification='All-cell exact search on integrated coordinates. Recall gates retain the existing HNSW qualification margins, fixed before this graph. No cell-type or treatment-preservation qualification.')
(a.out/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
