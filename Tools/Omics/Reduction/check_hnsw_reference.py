#!/usr/bin/env python3
"""Measure HNSW recall and independently verify every returned distance/fuzzy edge."""
import argparse,gzip,hashlib,importlib.metadata,json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
from umap.umap_ import fuzzy_simplicial_set
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--case',choices=['baron','hagai','norman'],required=True);p.add_argument('--protocol',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--exact-graph',type=Path);p.add_argument('--exact-scores',type=Path)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def data(path):return gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
protocol=json.loads(a.protocol.read_text());root=a.bundle;receipt=json.loads((root/'receipt.json').read_text())
for key,name in [('plan','plan.json'),('input','input/receipt.json'),('graph','graph.json'),('executionReport','execution.json')]:
    with (root/name).open('rb') as f:assert hashlib.file_digest(f,'sha256').digest()==bytes(receipt[key]['bytes']),name
plan=json.loads((root/'plan.json').read_text());execution=json.loads((root/'execution.json').read_text());h=execution['hnsw']
for key,value in protocol['parameters'].items():assert plan['approximation'][key]==value,(key,value)
assert h['options']==plan['approximation'] and execution['directedDistanceEvaluations']==h['constructionDistances']+h['queryDistances']
assert h['peakCachedScoreBytes']<=plan['approximation']['scoreCacheBytes']
g=json.loads((root/'graph.json').read_text());metadata=json.loads((root/'input/metadata.json').read_text());ir=json.loads((root/'input/receipt.json').read_text())
for key,name in [('metadata','metadata.json'),('scores','scores.bin')]:
    with (root/'input'/name).open('rb') as f:assert hashlib.file_digest(f,'sha256').digest()==bytes(ir[key]['bytes'])
assert g['method'].startswith('approximate-HNSW-') and g['distancePairs']==execution['directedDistanceEvaluations']
assert g['cells']==[dict(sampleID=c['sampleID'],barcode=c['barcode']) for c in metadata['cells']]
n=len(g['cells']);d=g['dimensions'];k=g['options']['neighbors'];assert n==dict(baron=8569,hagai=13863,norman=111445)[a.case] and k==15
path=root/'input/scores.bin';assert path.stat().st_size==n*d*16
records=np.fromfile(path,dtype=[('row','<u4'),('component','<u4'),('value','<f8')]);np.testing.assert_array_equal(records['row'],np.repeat(np.arange(n),d));np.testing.assert_array_equal(records['component'],np.tile(np.arange(d),n));x=records['value'].reshape(n,d).copy();assert np.all(np.isfinite(x));del records
indices=np.asarray(g['neighborIndices']).reshape(n,k);distances=np.asarray(g['neighborDistances']).reshape(n,k)
assert np.all((indices>=0)&(indices<n)) and np.all(np.isfinite(distances)) and np.all(distances>=0)
np.testing.assert_array_equal(indices[:,0],np.arange(n));assert np.all(distances[:,0]==0) and np.all(np.diff(np.sort(indices,axis=1),axis=1)!=0)
maxDistanceError=0.0
for start in range(0,n,1024):
    end=min(n,start+1024);delta=x[start:end,None,:]-x[indices[start:end]];expected=np.sqrt(np.sum(delta*delta,axis=2))
    np.testing.assert_allclose(distances[start:end],expected,rtol=1e-12,atol=1e-12)
    maxDistanceError=max(maxDistanceError,float(np.max(np.abs(distances[start:end]-expected))))
    for row in range(start,end):assert np.array_equal(np.lexsort((indices[row,1:],distances[row,1:])),np.arange(k-1))
if a.exact_graph:
    assert a.exact_scores is not None
    assert data(a.exact_scores)==path.read_bytes(),'PCA score bytes differ from exact-neighbor reference'
    exact=json.loads(data(a.exact_graph));assert exact['cells']==g['cells'] and exact['dimensions']==d
    assert exact['options']['neighbors']==k
    rows=np.arange(n);expectedIndices=np.asarray(exact['neighborIndices']).reshape(n,k)[:,1:];expectedDistances=np.asarray(exact['neighborDistances']).reshape(n,k)[:,1:]
    exactGraph=sparse.csr_matrix((exact['weights'],exact['columnIndices'],exact['rowOffsets']),shape=(n,n));del exact
else:
    assert a.case=='norman'
    rows=np.sort(np.random.default_rng(protocol['norman']['seed']).choice(n,protocol['norman']['uniformSampleCells'],replace=False))
    expectedIndices=np.empty((len(rows),k-1),dtype=np.int64);expectedDistances=np.empty((len(rows),k-1),dtype=np.float64)
    for start in range(0,len(rows),32):
        block=cdist(x[rows[start:start+32]],x,metric='euclidean')
        for slot,distance in enumerate(block):
            row=rows[start+slot];distance[row]=np.inf
            # Exact distance/index ordering, including all ties at the cutoff.
            cutoff=np.partition(distance,k-2)[k-2];candidates=np.flatnonzero(distance<=cutoff)
            chosen=candidates[np.lexsort((candidates,distance[candidates]))[:k-1]]
            expectedIndices[start+slot]=chosen;expectedDistances[start+slot]=distance[chosen]
    exactGraph=None
retrieved=indices[rows,1:]
recall=np.any(retrieved[:,:,None]==expectedIndices[:,None,:],axis=2).mean(axis=1)
# Report cutoff-distance equivalence separately; it never replaces the strict gate.
tieAware=np.mean(distances[rows,1:]<=expectedDistances[:,-1,None]+1e-12,axis=1)
np.savez_compressed(a.out/'recall-queries.npz',rows=rows,exactIndices=expectedIndices,exactDistances=expectedDistances,strictRecall=recall,tieAwareRecall=tieAware)
reference,sigma,rho=fuzzy_simplicial_set(x,n_neighbors=k,random_state=np.random.RandomState(7),metric='euclidean',knn_indices=indices,knn_dists=np.ascontiguousarray(distances,dtype=np.float32),local_connectivity=g['options']['localConnectivity'],set_op_mix_ratio=1.0)
reference=reference.tocsr();reference.eliminate_zeros();reference.sort_indices()
native=sparse.csr_matrix((g['weights'],g['columnIndices'],g['rowOffsets']),shape=(n,n))
assert (native-native.T).nnz==0 and np.all(native.diagonal()==0) and np.all((native.data>0)&(native.data<=1))
error=(native-reference).tocoo();maximumWeightError=float(np.max(np.abs(error.data))) if error.nnz else 0.0
referenceAtError=np.asarray(reference[error.row,error.col]).ravel()
assert np.all(np.abs(error.data)<=1e-5+1e-4*np.abs(referenceAtError))
np.testing.assert_allclose(g['rhos'],rho,rtol=1e-6,atol=1e-6);np.testing.assert_allclose(g['sigmas'],sigma,rtol=1e-4,atol=1e-6)
mass=np.exp(-np.maximum(distances[:,1:]-np.asarray(g['rhos'])[:,None],0)/np.asarray(g['sigmas'])[:,None]).sum(axis=1)
np.testing.assert_allclose(g['kernelMassResiduals'],np.abs(mass-np.log2(k)),rtol=1e-8,atol=1e-12)
components=sparse.csgraph.connected_components(native,directed=False,return_labels=False);assert components==g['connectedComponents']
assert int(np.count_nonzero(np.diff(native.indptr)==0))==g['isolatedCells']
gate=float(np.mean(recall))>=protocol['minimumMeanStrictRecall'] and float(np.quantile(recall,0.05))>=protocol['minimumFifthPercentileStrictRecall']
result=dict(status='passed' if gate else 'recall-gate-failed',case=a.case,cells=n,dimensions=d,queryCells=len(rows),allCellRecall=a.exact_graph is not None,meanStrictRecall=float(np.mean(recall)),fifthPercentileStrictRecall=float(np.quantile(recall,0.05)),minimumStrictRecall=float(np.min(recall)),meanTieAwareRecall=float(np.mean(tieAware)),maximumReturnedDistanceError=maxDistanceError,maximumFuzzyWeightError=maximumWeightError,fuzzyTopologyExact=bool(np.array_equal(native.indptr,reference.indptr) and np.array_equal(native.indices,reference.indices)),connectivityEntries=native.nnz,connectedComponents=components,distanceEvaluations=execution['directedDistanceEvaluations'],peakCachedScoreBytes=h['peakCachedScoreBytes'],indexStorageBytes=h['indexStorageBytes'],protocolSHA256=hashlib.sha256(a.protocol.read_bytes()).hexdigest(),graphSHA256=hashlib.sha256((root/'graph.json').read_bytes()).hexdigest(),versions={v:importlib.metadata.version(v) for v in ['numpy','scipy','umap-learn']},scope='Numerical returned-distance/fuzzy-graph and declared exact-query recall checks; no biological generalization or million-cell qualification')
if exactGraph is not None:
    support=native.copy();support.data[:]=1;other=exactGraph.copy();other.data[:]=1
    intersection=support.multiply(other).nnz;result['exactFuzzyEdgeJaccard']=intersection/(support.nnz+other.nnz-intersection)
(a.out/'checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result),flush=True)
assert gate,'Predeclared strict recall gate failed; retain this result before changing parameters'
