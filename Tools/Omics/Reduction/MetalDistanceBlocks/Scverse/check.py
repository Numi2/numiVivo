from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import numpy as np,json,hashlib,statistics,importlib.metadata as md
from scipy import sparse
r=Path(__file__).parent;p=json.loads((r/'protocol.json').read_text());term=json.loads((r/'terminal.json').read_text());assert len(term['results'])==9
native={};scverse=None
for entry in term['results']:
 folder=r/(str(entry['run'])+'-'+entry['backend'])
 if entry['backend']=='scanpy':
  matrices=[sparse.load_npz(folder/name).tocsr() for name in ['distances.npz','connectivities.npz']]
  for m in matrices:m.sort_indices();assert m.shape==(24673,24673) and np.all(np.isfinite(m.data))
  if scverse is None:scverse=matrices
  else:
   for a,b in zip(scverse,matrices):assert np.array_equal(a.indptr,b.indptr) and np.array_equal(a.indices,b.indices) and np.array_equal(a.data,b.data)
 else:
  g=json.loads((folder/'graph.json').read_text());assert len(g['cells'])==24673
  if entry['backend'] in native:assert native[entry['backend']]==g
  else:native[entry['backend']]=g
assert native['cpu']['cells']==native['metal']['cells']
comparisons={}
for mode,g in native.items():
 ids=np.asarray(g['neighborIndices']).reshape(24673,20);values=np.asarray(g['neighborDistances']).reshape(24673,20)
 missing=0;error=0.;referenceEntries=0
 dist,conn=scverse
 for row in range(24673):
  lo,hi=dist.indptr[row:row+2];js=dist.indices[lo:hi];ds=dist.data[lo:hi];reference={int(j):float(v) for j,v in zip(js,ds) if j!=row}
  candidate={int(j):float(v) for j,v in zip(ids[row,1:],values[row,1:])};referenceEntries+=len(reference);missing+=len(set(reference)-set(candidate))
  if set(reference)&set(candidate):error=max(error,max(abs(reference[j]-candidate[j]) for j in set(reference)&set(candidate)))
 graph=sparse.csr_matrix((g['weights'],g['columnIndices'],g['rowOffsets']),shape=conn.shape);graph.sort_indices();keys=np.repeat(np.arange(24673,dtype=np.int64),np.diff(graph.indptr))*24673+graph.indices;sk=np.repeat(np.arange(24673,dtype=np.int64),np.diff(conn.indptr))*24673+conn.indices
 common,ig,isc=np.intersect1d(keys,sk,assume_unique=True,return_indices=True)
 comparisons[mode]=dict(scanpyNonSelfNeighborEntries=referenceEntries,missingScanpyNeighbors=missing,neighborRecall=1-missing/referenceEntries,maximumCommonNeighborDistanceDifference=error,nativeEdges=len(keys),scanpyEdges=len(sk),commonEdges=len(common),maximumCommonFuzzyWeightDifference=float(abs(graph.data[ig]-conn.data[isc]).max()))
summary={mode:dict(graphStageSeconds=[x['graphStageSeconds'] for x in term['results'] if x['backend']==mode],medianGraphStageSeconds=statistics.median(x['graphStageSeconds'] for x in term['results'] if x['backend']==mode),processWallSeconds=[x['processWallSeconds'] for x in term['results'] if x['backend']==mode],peakRSSBytes=[x['peakRSSBytes'] for x in term['results'] if x['backend']==mode]) for mode in ['cpu','metal','scanpy']}
result=dict(scope=p['scope'],cells=24673,repetitionsExact=True,comparisons=comparisons,timings=summary,packages={n:md.version(n) for n in ['scanpy','umap-learn','numpy','scipy','scikit-learn','numba','anndata']},protocolSHA256=hashlib.sha256((r/'protocol.json').read_bytes()).hexdigest(),fullCLISpeedupClaim=False,biologicalQualification=False)
(r/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
