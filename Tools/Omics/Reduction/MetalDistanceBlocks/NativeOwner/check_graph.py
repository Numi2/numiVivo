from pathlib import Path
import numpy as np,json
r=Path(__file__).parent;record=np.dtype([('row','<u4'),('col','<u4'),('value','<f8')]);n=13863;k=20
cpu=np.fromfile(r/'cpu/neighbors.bin',record);metal=np.fromfile(r/'metal/neighbors.bin',record);oracle=np.fromfile('/Users/n/numivivo-metal-knn-hoisted-20260912/1-metal/indices.bin','<u4')
assert np.array_equal(metal['col'],oracle)
a=cpu['col'].reshape(n,k);b=metal['col'].reshape(n,k);assert all(set(x)==set(y) for x,y in zip(a,b))
reference=json.loads((r/'metal-json/graph.json').read_text());assert reference['method']=='metal-FP32-windowed-PCA-knn-v1'
assert np.array_equal(metal['col'],reference['neighborIndices']);np.testing.assert_array_equal(metal['value'],reference['neighborDistances'])
edges=np.fromfile(r/'metal/edges.bin',record);cpu_edges=np.fromfile(r/'cpu/edges.bin',record)
assert np.array_equal(edges['row'],cpu_edges['row']) and np.array_equal(edges['col'],cpu_edges['col'])
assert np.array_equal(edges['col'],reference['columnIndices'])
assert np.array_equal(edges['row'],np.repeat(np.arange(n),np.diff(reference['rowOffsets'])))
np.testing.assert_allclose(edges['value'],reference['weights'],atol=1e-12,rtol=1e-12)
assert np.all(np.isfinite(edges['value'])) and np.all((edges['value']>0)&(edges['value']<=1))
weightsError=float(abs(edges['value']-cpu_edges['value']).max())
# Precision qualification is measured, not silently promoted to biology.
reports={mode:json.loads((r/mode/'graph.json').read_text()) for mode in ['cpu','metal']}
assert reports['cpu']['connectedComponents']==reports['metal']['connectedComponents']
assert reports['cpu']['isolatedCells']==reports['metal']['isolatedCells']
result=dict(cells=n,neighborEntries=len(metal),allNeighborMembershipExact=True,rowsWithOrderDifference=int(np.count_nonzero(np.any(a!=b,axis=1))),metalBinaryJSONNeighborsExact=True,connectivityEntries=len(edges),allEdgeCoordinatesExact=True,maximumFuzzyWeightDifferenceVsCPU=weightsError,maximumBinaryJSONWeightDifference=float(abs(edges['value']-reference['weights']).max()),connectedComponents=reports['metal']['connectedComponents'],isolatedCells=reports['metal']['isolatedCells'],biologicalQualification=False)
(r/'graph-check.json').write_text(json.dumps(result,indent=2));print(result)
