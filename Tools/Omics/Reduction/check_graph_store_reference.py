#!/usr/bin/env python3
"""Every binary record versus independently qualified JSON graph and PCA scores."""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
from scipy import sparse
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--baseline-graph',type=Path,required=True);p.add_argument('--baseline-scores',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def raw(path):return gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').digest()
r=a.bundle;receipt=json.loads((r/'receipt.json').read_text());g=json.loads((r/'graph.json').read_text());plan=json.loads((r/'plan.json').read_text());base=json.loads(raw(a.baseline_graph))
assert plan['storage']=='binary'
for key,path in [('plan','plan.json'),('graph','graph.json'),('executionReport','execution.json'),('input','input/receipt.json')]:assert sha(r/path)==bytes(receipt[key]['bytes'])
for key,name in [('neighbors','neighbors.bin'),('bandwidths','bandwidths.bin'),('offsets','offsets.bin'),('edges','edges.bin')]:assert sha(r/name)==bytes(g[key]['bytes'])
ir=json.loads((r/'input/receipt.json').read_text());assert sha(r/'input/scores.bin')==bytes(ir['scores']['bytes']);assert sha(r/'input/metadata.json')==bytes(ir['metadata']['bytes'])
assert (r/'input/scores.bin').read_bytes()==raw(a.baseline_scores)
metadata=json.loads((r/'input/metadata.json').read_text());assert base['cells']==[dict(sampleID=c['sampleID'],barcode=c['barcode']) for c in metadata['cells']]
n=g['cells'];k=g['options']['neighbors'];assert n==len(base['cells']) and g['dimensions']==base['dimensions'] and g['method']==base['method']
for key in ['options','connectedComponents','isolatedCells','distancePairs']:assert g[key]==base[key],key
record=np.dtype([('row','<u4'),('column','<u4'),('bits','<u8')])
def records(name,count):
    path=r/name;assert path.stat().st_size==count*16
    return np.fromfile(path,dtype=record)
v=records('neighbors.bin',n*k)
np.testing.assert_array_equal(v['row'],np.repeat(np.arange(n),k));np.testing.assert_array_equal(v['column'],base['neighborIndices']);np.testing.assert_array_equal(v['bits'],np.asarray(base['neighborDistances'],dtype='<f8').view('<u8'))
v=records('bandwidths.bin',n*3)
np.testing.assert_array_equal(v['row'],np.repeat(np.arange(n),3));np.testing.assert_array_equal(v['column'],np.tile(np.arange(3),n))
np.testing.assert_array_equal(v['bits'],np.column_stack([base['rhos'],base['sigmas'],base['kernelMassResiduals']]).astype('<f8').ravel().view('<u8'))
v=records('offsets.bin',n+1)
np.testing.assert_array_equal(v['row'],np.arange(n+1));assert np.all(v['column']==0);np.testing.assert_array_equal(v['bits'],base['rowOffsets']);offsets=v['bits'].astype(np.int64)
v=records('edges.bin',g['connectivityEntries'])
np.testing.assert_array_equal(v['row'],np.repeat(np.arange(n),np.diff(offsets)));np.testing.assert_array_equal(v['column'],base['columnIndices']);np.testing.assert_array_equal(v['bits'],np.asarray(base['weights'],dtype='<f8').view('<u8'))
csr=sparse.csr_matrix((v['bits'].view('<f8'),v['column'],offsets),shape=(n,n));assert (csr-csr.T).nnz==0 and not csr.diagonal().any();assert sparse.csgraph.connected_components(csr,directed=False,return_labels=False)==g['connectedComponents']
assert not any(p.name.startswith('.graph-transpose-') for p in r.iterdir())
result=dict(status='passed',cells=n,dimensions=g['dimensions'],neighborEntries=n*k,connectivityEntries=g['connectivityEntries'],allBinaryCoordinatesAndFP64BitsExact=True,allInputScoreBytesExact=True,allCellIdentitiesExact=True,components=g['connectedComponents'],baselineGraphSHA256=hashlib.sha256(raw(a.baseline_graph)).hexdigest(),baselineScoresSHA256=hashlib.sha256(raw(a.baseline_scores)).hexdigest(),graphSHA256=sha(r/'graph.json').hex(),binaryBytes=sum((r/name).stat().st_size for name in ['neighbors.bin','bandwidths.bin','offsets.bin','edges.bin']),scope='Every record equals previously independently qualified full graph; preserves its recall scope, not new biological or million-cell qualification')
(a.out/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
