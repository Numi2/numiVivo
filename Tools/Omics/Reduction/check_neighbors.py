#!/usr/bin/env python3
"""Independent chunked exact kNN and umap-learn fuzzy graph comparison."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
from umap.umap_ import fuzzy_simplicial_set
p = argparse.ArgumentParser()
p.add_argument('--report', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
r = json.loads(a.report.read_text())
g, reduction = r['neighbors'], r['reduction']
assert g['cells'] == reduction['cells']
x = np.asarray(reduction['scores'], dtype=np.float64)
n, dimensions = x.shape
assert g['dimensions'] == dimensions
k = g['options']['neighbors']
indices = np.asarray(g['neighborIndices']).reshape(n,k)
distances = np.asarray(g['neighborDistances']).reshape(n,k)
assert np.array_equal(indices[:,0], np.arange(n))
assert np.all(distances[:,0] == 0)
maximum_distance_error = 0.0
# At most 64-by-cells distances, never a dense cells-by-cells matrix.
for start in range(0,n,64):
    block = cdist(x[start:start+64], x, metric='euclidean')
    for slot, row in enumerate(block):
        i = start+slot
        order = np.lexsort((np.arange(n),row))
        expected = order[order != i][:k-1]
        assert np.array_equal(indices[i,1:], expected), ('neighbor order',i)
        np.testing.assert_allclose(distances[i,1:],row[expected],rtol=1e-12,atol=1e-12)
        maximum_distance_error = max(maximum_distance_error,float(np.max(np.abs(distances[i,1:]-row[expected]))))
reference, sigma, rho = fuzzy_simplicial_set(x,n_neighbors=k,random_state=np.random.RandomState(7),metric='euclidean',
    knn_indices=indices,knn_dists=np.ascontiguousarray(distances,dtype=np.float32),
    local_connectivity=g['options']['localConnectivity'],set_op_mix_ratio=1.0)
reference = reference.tocsr(); reference.eliminate_zeros(); reference.sort_indices()
native = sparse.csr_matrix((g['weights'],g['columnIndices'],g['rowOffsets']),shape=(n,n))
assert np.array_equal(native.indptr,reference.indptr)
assert np.array_equal(native.indices,reference.indices)
np.testing.assert_allclose(native.data,reference.data,rtol=1e-4,atol=1e-5)
np.testing.assert_allclose(g['rhos'],rho,rtol=1e-6,atol=1e-6)
np.testing.assert_allclose(g['sigmas'],sigma,rtol=1e-4,atol=1e-6)
delta = distances[:,1:]-np.asarray(g['rhos'])[:,None]
mass = np.exp(-np.maximum(delta,0)/np.asarray(g['sigmas'])[:,None]).sum(axis=1)
np.testing.assert_allclose(g['kernelMassResiduals'],np.abs(mass-np.log2(k)),atol=1e-12,rtol=1e-8)
assert (native-native.T).nnz == 0
assert np.all(native.diagonal() == 0)
assert np.all((native.data>0)&(native.data<=1))
components = sparse.csgraph.connected_components(native,directed=False,return_labels=False)
assert components == g['connectedComponents']
assert int(np.count_nonzero(np.diff(native.indptr)==0)) == g['isolatedCells']
assert g['distancePairs'] == n*(n-1)//2
result = dict(status='passed',reportSHA256=hashlib.sha256(a.report.read_bytes()).hexdigest(),cells=n,dimensions=dimensions,
    neighborsIncludingSelf=k,distancePairs=g['distancePairs'],connectivityEntries=native.nnz,connectedComponents=components,
    exactNeighborMembershipAndOrder=True,maximumDistanceError=maximum_distance_error,
    maximumConnectivityError=float(np.max(np.abs(native.data-reference.data))),
    maximumSigmaError=float(np.max(np.abs(np.asarray(g['sigmas'])-sigma))),
    maximumKernelMassResidual=max(g['kernelMassResiduals']),
    qualification='Numerical exact-neighbor and fuzzy-graph agreement; no embedding, clustering or biological validation',
    versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','umap-learn']})
a.out.parent.mkdir(parents=True,exist_ok=True)
a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(result,indent=2))
