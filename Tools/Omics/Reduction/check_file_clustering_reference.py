#!/usr/bin/env python3
"""Full sparse objective/connectivity, frozen oracle, and three Louvain references."""
import argparse,hashlib,importlib.metadata,json
from pathlib import Path
import numpy as np
import networkx as nx
from scipy import sparse
from sklearn.metrics import adjusted_rand_score,adjusted_mutual_info_score,homogeneity_completeness_v_measure
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--oracle-result',type=Path,help='Optional frozen result; the native runner owns this comparison when omitted');p.add_argument('--protocol',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').digest()
r=a.bundle;protocol=json.loads(a.protocol.read_text());receipt=json.loads((r/'receipt.json').read_text())
for key,name in [('input','input/receipt.json'),('plan','plan.json'),('result','result.json'),('executionReport','execution.json')]:assert sha(r/name)==bytes(receipt[key]['bytes'])
if a.oracle_result is not None:assert (r/'result.json').read_bytes()==a.oracle_result.read_bytes(),'Frozen pre-refactor result differs'
c=json.loads((r/'result.json').read_text());plan=json.loads((r/'plan.json').read_text());execution=json.loads((r/'execution.json').read_text());root=r/'input';g=json.loads((root/'graph.json').read_text());gr=json.loads((root/'receipt.json').read_text());ir=json.loads((root/'input/receipt.json').read_text())
for key,name in [('input','input/receipt.json'),('plan','plan.json'),('graph','graph.json'),('executionReport','execution.json')]:assert sha(root/name)==bytes(gr[key]['bytes'])
for key in ['neighbors','bandwidths','offsets','edges']:assert sha(root/(key+'.bin'))==bytes(g[key]['bytes'])
assert sha(root/'input/metadata.json')==bytes(ir['metadata']['bytes'])
metadata=json.loads((root/'input/metadata.json').read_text());assert c['cells']==[dict(sampleID=v['sampleID'],barcode=v['barcode']) for v in metadata['cells']]
assert c['options']==plan['clustering'] and c['options']['resolution']==protocol['resolution'] and c['options']['seed']==protocol['nativeSeed']
assert 0<execution['edgeVisits']<=plan['maximumEdgeVisits'] and execution['aggregatedLevels']==len(c['levels'])-1
n=len(c['cells']);assert n==g['cells'];labels=np.asarray(c['labels'],dtype=np.int64);assert labels.shape==(n,)
record=np.dtype([('row','<u4'),('column','<u4'),('bits','<u8')]);o=np.fromfile(root/'offsets.bin',dtype=record);e=np.fromfile(root/'edges.bin',dtype=record)
assert len(o)==n+1 and len(e)==g['connectivityEntries'];np.testing.assert_array_equal(o['row'],np.arange(n+1));assert np.all(o['column']==0)
offsets=o['bits'].astype(np.int64);assert offsets[0]==0 and offsets[-1]==len(e) and np.all(np.diff(offsets)>=0)
np.testing.assert_array_equal(e['row'],np.repeat(np.arange(n),np.diff(offsets)));w=sparse.csr_matrix((e['bits'].view('<f8'),e['column'],offsets),shape=(n,n));assert (w-w.T).nnz==0 and not w.diagonal().any() and np.all(w.data>0)
graph=nx.from_scipy_sparse_array(w)
partition=[set(np.flatnonzero(labels==j).tolist()) for j in range(len(c['clusterSizes']))]
assert [len(v) for v in partition]==c['clusterSizes'] and np.unique(labels).tolist()==list(range(len(partition)))
assert [min(v) for v in partition]==sorted(min(v) for v in partition)
objective=nx.community.modularity(graph,partition,weight='weight',resolution=protocol['resolution']);np.testing.assert_allclose(c['modularity'],objective,rtol=1e-10,atol=1e-10)
disconnected=sum(not nx.is_connected(graph.subgraph(v)) for v in partition);assert disconnected==c['disconnectedCommunities']
assert all(b['modularity']+1e-10>=a['modularity'] for a,b in zip(c['levels'],c['levels'][1:]))
assert c['levels'][-1]['communities']==len(partition)
references=[];arrays=dict(native=labels)
for seed in protocol['referenceSeeds']:
    groups=nx.community.louvain_communities(graph,weight='weight',resolution=protocol['resolution'],threshold=c['options']['levelTolerance'],seed=seed)
    values=np.empty(n,dtype=np.int64)
    for j,group in enumerate(groups):values[list(group)]=j
    arrays['seed'+str(seed)]=values
    references.append(dict(seed=seed,clusters=len(groups),modularity=nx.community.modularity(graph,groups,weight='weight',resolution=protocol['resolution']),adjustedRand=float(adjusted_rand_score(labels,values))))
    (a.out/'references-progress.json').write_text(json.dumps(references,indent=2)+'\n');print('reference seed '+str(seed)+' completed',flush=True)
np.savez_compressed(a.out/'partitions.npz',**arrays)
samples={v['id']:v for v in metadata['samples']};associations={}
for name,values in [('annotatedCellGroup',[v.get('group') for v in metadata['cells']]),('donor',[samples[v['sampleID']].get('donorID') for v in metadata['cells']]),('condition',[samples[v['sampleID']].get('condition') for v in metadata['cells']])]:
    complete=all(v is not None for v in values);available=complete and len(set(values))>1
    entry=dict(status='descriptive' if available else 'missing-or-single-level',complete=complete,levels=len(set(values)))
    if available:
        h,cpl,v=homogeneity_completeness_v_measure(values,labels)
        entry.update(adjustedRand=float(adjusted_rand_score(values,labels)),adjustedMutualInformation=float(adjusted_mutual_info_score(values,labels)),homogeneity=h,completeness=cpl,vMeasure=v)
    associations[name]=entry
nativePairwise=[float(adjusted_rand_score(arrays['seed'+str(x)],arrays['seed'+str(y)])) for i,x in enumerate(protocol['referenceSeeds']) for y in protocol['referenceSeeds'][i+1:]]
gap=max(v['modularity'] for v in references)-objective;gate=gap<=protocol['maximumReferenceModularityDeficit']
result=dict(status='passed' if gate else 'reference-quality-gap',cells=n,clusters=len(partition),allFrozenLegacyResultBytesExact=True if a.oracle_result is not None else None,nativeModularity=c['modularity'],independentModularity=objective,disconnectedCommunities=disconnected,referenceMaximumModularityGap=gap,references=references,referencePairwiseARI=nativePairwise,metadataAssociations=associations,resultSHA256=sha(r/'result.json').hex(),graphSHA256=sha(root/'graph.json').hex(),protocolSHA256=sha(a.protocol).hex(),edgeVisits=execution['edgeVisits'],versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','networkx','scikit-learn']},scope='Full graph numerical objective, connectivity and heuristic comparison; descriptive metadata association is not authoritative cell annotation, donor integration or biological qualification')
(a.out/'checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result),flush=True)
assert gate,'Predeclared reference modularity deficit exceeded; retain and investigate'
