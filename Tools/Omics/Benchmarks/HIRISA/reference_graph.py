#!/usr/bin/env python3
"""Frozen exact query panel, then full binary HNSW distance/fuzzy-graph checks."""
import argparse
import hashlib
import importlib.metadata
import json
import time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist

RECORD = np.dtype([('row', '<u4'), ('column', '<u4'), ('value', '<f8')])

def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')

def matrix(path, rows, columns):
    assert path.stat().st_size == rows * columns * 16
    records = np.memmap(path, mode='r', dtype=RECORD, shape=(rows, columns))
    for start in range(0, rows, 8192):
        block = records[start:start + 8192]
        np.testing.assert_array_equal(block['row'], np.broadcast_to(np.arange(start, start + len(block))[:, None], block.shape))
        np.testing.assert_array_equal(block['column'], np.broadcast_to(np.arange(columns), block.shape))
        assert np.isfinite(block['value']).all()
    return records['value'].copy()

def prepare(a, protocol):
    start = time.monotonic()
    n, d, k = (protocol[v] for v in ['cells', 'dimensions', 'neighbors'])
    assert sha(a.scores) == protocol['scoresSHA256']
    x = matrix(a.scores, n, d)
    rows = np.sort(np.random.default_rng(protocol['querySeed']).choice(n, protocol['queryCells'], replace=False))
    indices = np.empty((len(rows), k-1), dtype=np.uint32)
    distances = np.empty((len(rows), k-1), dtype=np.float64)
    # Only 16 query-by-cell distances are resident at once. Every cell is a
    # candidate; exact cutoff ties are resolved by original row identity.
    for first in range(0, len(rows), 16):
        block = cdist(x[rows[first:first+16]], x, metric='euclidean')
        for slot, distance in enumerate(block):
            distance[rows[first+slot]] = np.inf
            cutoff = np.partition(distance, k-2)[k-2]
            candidates = np.flatnonzero(distance <= cutoff)
            chosen = candidates[np.lexsort((candidates, distance[candidates]))[:k-1]]
            indices[first+slot] = chosen
            distances[first+slot] = distance[chosen]
        if first % 256 == 0:
            print('exact queries', first + len(block), '/', len(rows), flush=True)
    np.savez_compressed(a.out/'exact-queries.npz', rows=rows, indices=indices, distances=distances)
    write(a.out/'reference.json', dict(status='passed', protocolSHA256=sha(a.protocol), scoresSHA256=sha(a.scores),
          queriesSHA256=sha(a.out/'exact-queries.npz'), cells=n, queryCells=len(rows),
          candidateDistancePairs=n*len(rows), maximumDistanceBlockBytes=16*n*8,
          seconds=time.monotonic()-start, versions={v:importlib.metadata.version(v) for v in ['numpy','scipy']},
          scope='Exact FP64 SciPy distances for frozen queries against every original PCA row; not all-cell exact recall.'))

def check(a, protocol):
    from umap.umap_ import fuzzy_simplicial_set
    start = time.monotonic()
    root = a.bundle
    n, d, k = (protocol[v] for v in ['cells', 'dimensions', 'neighbors'])
    receipt = json.loads((root/'receipt.json').read_text())
    for key, path in [('plan','plan.json'), ('input','input/receipt.json'), ('graph','graph.json'), ('executionReport','execution.json')]:
        assert sha(root/path) == bytes(receipt[key]['bytes']).hex(), path
    plan = json.loads((root/'plan.json').read_text())
    # Codable fills omitted defaults. Compare the frozen explicit settings.
    for key in ['schemaVersion','inputKind','storage']:
        assert plan[key] == protocol['plan'][key]
    for group in ['neighbors','execution','approximation']:
        for key, value in protocol['plan'][group].items():
            assert plan[group][key] == value, (group, key)
    graph = json.loads((root/'graph.json').read_text())
    execution = json.loads((root/'execution.json').read_text())
    assert execution['cells'] == n and execution['dimensions'] == d
    assert graph['cells'] == n and graph['dimensions'] == d and graph['neighborEntries'] == n*k
    assert graph['method'].startswith('approximate-HNSW-') and graph['options'] == plan['neighbors']
    hnsw = execution['hnsw']
    assert hnsw['options'] == plan['approximation']
    assert graph['distancePairs'] == execution['directedDistanceEvaluations'] == hnsw['constructionDistances'] + hnsw['queryDistances']
    assert execution['scalarDistanceTerms'] == graph['distancePairs'] * d
    assert graph['distancePairs'] <= plan['approximation']['maximumDistanceEvaluations']
    assert hnsw['peakCachedScoreBytes'] <= plan['approximation']['scoreCacheBytes']
    input_receipt = json.loads((root/'input/receipt.json').read_text())
    assert bytes(input_receipt['source']['bytes']).hex() == protocol['sourceSHA256']
    assert input_receipt['implementation'] == receipt['implementation']
    for name in ['scores','metadata','quality']:
        path = root/'input'/(name + ('.bin' if name == 'scores' else '.json'))
        assert sha(path) == bytes(input_receipt[name]['bytes']).hex() == protocol[name+'SHA256']
    for key in ['neighbors','bandwidths','offsets','edges']:
        assert sha(root/(key+'.bin')) == bytes(graph[key]['bytes']).hex(), key
    reference = json.loads((a.reference/'reference.json').read_text())
    assert reference['status'] == 'passed' and reference['protocolSHA256'] == sha(a.protocol)
    assert reference['scoresSHA256'] == protocol['scoresSHA256']
    assert reference['queriesSHA256'] == sha(a.reference/'exact-queries.npz')
    queries = np.load(a.reference/'exact-queries.npz', allow_pickle=False)
    np.testing.assert_array_equal(queries['rows'], np.sort(np.random.default_rng(protocol['querySeed']).choice(n, protocol['queryCells'], replace=False)))
    x = matrix(root/'input/scores.bin', n, d)
    assert (root/'neighbors.bin').stat().st_size == n*k*16
    records = np.memmap(root/'neighbors.bin', mode='r', dtype=RECORD, shape=(n,k))
    indices, distances = records['column'], records['value']
    maximum_distance_error = 0.0
    for first in range(0,n,1024):
        end = min(n,first+1024)
        index, distance = indices[first:end], distances[first:end]
        np.testing.assert_array_equal(records['row'][first:end], np.broadcast_to(np.arange(first,end)[:,None], index.shape))
        assert (index<n).all() and np.isfinite(distance).all() and (distance>=0).all()
        np.testing.assert_array_equal(index[:,0], np.arange(first,end))
        assert (distance[:,0]==0).all() and (np.diff(np.sort(index,axis=1),axis=1)>0).all()
        delta = x[first:end,None,:]-x[index]
        expected = np.sqrt(np.sum(delta*delta,axis=2))
        np.testing.assert_allclose(distance,expected,rtol=protocol['distanceRelativeTolerance'],atol=protocol['distanceAbsoluteTolerance'])
        maximum_distance_error = max(maximum_distance_error,float(np.max(np.abs(distance-expected))))
        ordered = (distance[:,2:]>distance[:,1:-1]) | ((distance[:,2:]==distance[:,1:-1]) & (index[:,2:]>index[:,1:-1]))
        assert ordered.all()
    retrieved = indices[queries['rows'],1:]
    recall = np.any(retrieved[:,:,None]==queries['indices'][:,None,:],axis=2).mean(axis=1)
    tie_aware = np.mean(distances[queries['rows'],1:]<=queries['distances'][:,-1,None]+1e-12,axis=1)
    np.savez_compressed(a.out/'recall.npz',rows=queries['rows'],strictRecall=recall,tieAwareRecall=tie_aware)
    print('All returned distances passed; mean strict recall',float(recall.mean()),flush=True)
    # Full umap-learn reference on the returned neighbors, not a second ANN fit.
    fuzzy,sigma,rho = fuzzy_simplicial_set(x,n_neighbors=k,random_state=np.random.RandomState(7),metric='euclidean',
        knn_indices=np.ascontiguousarray(indices,dtype=np.int32),knn_dists=np.ascontiguousarray(distances,dtype=np.float32),
        local_connectivity=plan['neighbors']['localConnectivity'],set_op_mix_ratio=1.0)
    fuzzy = fuzzy.tocsr(); fuzzy.eliminate_zeros(); fuzzy.sort_indices()
    bandwidths = matrix(root/'bandwidths.bin',n,3)
    np.testing.assert_allclose(bandwidths[:,0],rho,rtol=1e-6,atol=1e-6)
    np.testing.assert_allclose(bandwidths[:,1],sigma,rtol=1e-4,atol=1e-6)
    for first in range(0,n,8192):
        b=bandwidths[first:first+8192]
        mass=np.exp(-np.maximum(distances[first:first+8192,1:]-b[:,0,None],0)/b[:,1,None]).sum(axis=1)
        np.testing.assert_allclose(b[:,2],np.abs(mass-np.log2(k)),rtol=1e-8,atol=1e-12)
    integer_record=np.dtype([('row','<u4'),('column','<u4'),('bits','<u8')])
    assert (root/'offsets.bin').stat().st_size == (n+1)*16
    offsets_file=np.memmap(root/'offsets.bin',mode='r',dtype=integer_record)
    np.testing.assert_array_equal(offsets_file['row'],np.arange(n+1));assert (offsets_file['column']==0).all()
    offsets=offsets_file['bits'].astype(np.int64)
    entries=graph['connectivityEntries']
    assert offsets[0]==0 and offsets[-1]==entries and (np.diff(offsets)>=0).all()
    assert (root/'edges.bin').stat().st_size==entries*16
    edges=np.memmap(root/'edges.bin',mode='r',dtype=RECORD)
    for first in range(0,n,8192):
        end=min(n,first+8192);block=edges[offsets[first]:offsets[end]]
        np.testing.assert_array_equal(block['row'],np.repeat(np.arange(first,end),np.diff(offsets[first:end+1])))
        assert (block['column']<n).all() and (block['column']!=block['row']).all()
        assert np.isfinite(block['value']).all() and ((block['value']>0)&(block['value']<=1)).all()
        keys=block['row'].astype(np.uint64)*n+block['column']
        assert (keys[1:]>keys[:-1]).all()
    native=sparse.csr_matrix((edges['value'],edges['column'],offsets),shape=(n,n))
    assert (native-native.T).nnz==0
    difference=(native-fuzzy).tocoo()
    at_error=np.asarray(fuzzy[difference.row,difference.col]).ravel()
    assert np.all(np.abs(difference.data)<=protocol['fuzzyAbsoluteTolerance']+protocol['fuzzyRelativeTolerance']*np.abs(at_error))
    components=sparse.csgraph.connected_components(native,directed=False,return_labels=False)
    assert components==graph['connectedComponents']
    assert np.count_nonzero(np.diff(offsets)==0)==graph['isolatedCells']
    gate=float(recall.mean())>=protocol['minimumMeanStrictRecall'] and float(np.quantile(recall,.05))>=protocol['minimumFifthPercentileStrictRecall']
    result=dict(status='passed' if gate else 'recall-gate-failed',cells=n,dimensions=d,neighbors=k,
        queryCells=len(recall),allCellRecall=False,meanStrictRecall=float(recall.mean()),fifthPercentileStrictRecall=float(np.quantile(recall,.05)),minimumStrictRecall=float(recall.min()),meanTieAwareRecall=float(tie_aware.mean()),
        allReturnedDistancesChecked=True,allBinaryCoordinatesChecked=True,allSourceScoreMetadataQualityBytesExact=True,
        maximumReturnedDistanceError=maximum_distance_error,maximumFuzzyWeightError=float(np.max(np.abs(difference.data))) if difference.nnz else 0.0,
        fuzzyTopologyExact=bool(np.array_equal(native.indptr,fuzzy.indptr) and np.array_equal(native.indices,fuzzy.indices)),
        connectivityEntries=entries,connectedComponents=components,isolatedCells=graph['isolatedCells'],
        seconds=time.monotonic()-start,protocolSHA256=sha(a.protocol),receiptSHA256=sha(root/'receipt.json'),
        versions={v:importlib.metadata.version(v) for v in ['numpy','scipy','umap-learn','numba']},scope=protocol['scope'])
    write(a.out/'checks.json',result);print(json.dumps(result),flush=True)
    assert gate,'Frozen recall gate failed; retain this result before any parameter change'

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['prepare','check']);p.add_argument('--protocol',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--scores',type=Path)
    p.add_argument('--bundle',type=Path);p.add_argument('--reference',type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    protocol=json.loads(a.protocol.read_text())
    if a.mode=='prepare':
        assert a.scores is not None
        prepare(a,protocol)
    else:
        assert a.bundle is not None and a.reference is not None
        check(a,protocol)
