#!/usr/bin/env python3
"""Independent modularity, connectivity and NetworkX Louvain comparisons."""
import argparse,hashlib,importlib.metadata,json
from pathlib import Path
import numpy as np
import networkx as nx
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
p=argparse.ArgumentParser()
p.add_argument('--report',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
r=json.loads(a.report.read_text()); c=r['clustering']; g=r['neighbors']
labels=np.asarray(c['labels']); n=len(labels); gamma=c['options']['resolution']
assert c['cells']==g['cells'] and n==len(g['cells'])
w=sparse.csr_matrix((g['weights'],g['columnIndices'],g['rowOffsets']),shape=(n,n))
graph=nx.from_scipy_sparse_array(w)
partition=[set(np.flatnonzero(labels==j).tolist()) for j in range(len(c['clusterSizes']))]
assert [len(s) for s in partition]==c['clusterSizes']
assert np.unique(labels).tolist()==list(range(len(partition)))
assert [min(s) for s in partition]==sorted(min(s) for s in partition)
objective=nx.community.modularity(graph,partition,weight='weight',resolution=gamma)
np.testing.assert_allclose(c['modularity'],objective,atol=1e-10,rtol=1e-10)
disconnected=sum(not nx.is_connected(graph.subgraph(s)) for s in partition)
assert disconnected==c['disconnectedCommunities']
assert all(b['modularity']+1e-10>=a['modularity'] for a,b in zip(c['levels'],c['levels'][1:]))
assert c['levels'][-1]['communities']==len(partition)
references=[]
for seed in [7,19,41]:
 groups=nx.community.louvain_communities(graph,weight='weight',resolution=gamma,threshold=c['options']['levelTolerance'],seed=seed)
 reference=np.zeros(n,dtype=int)
 for j,group in enumerate(groups): reference[list(group)]=j
 references.append(dict(seed=seed,clusters=len(groups),modularity=nx.community.modularity(graph,groups,weight='weight',resolution=gamma),
                        adjustedRand=float(adjusted_rand_score(labels,reference))))
# Partition identity is not an oracle: Louvain local optima depend on ordering.
# Retain every comparison and flag a material objective deficit rather than
# asserting identical labels or treating either result as true cell types.
gap=max(item['modularity'] for item in references)-objective
samples={sample['id']:sample for sample in r['processed']['dataset']['samples']}
conditions=[samples[cell['sampleID']].get('condition') for cell in c['cells']]
donors=[samples[cell['sampleID']].get('donorID') for cell in c['cells']]
associations={}
for name, values in [('condition',conditions),('donor',donors)]:
    available=all(value is not None for value in values) and len(set(values))>1
    associations[name]=dict(status='descriptive' if available else 'not-applicable-missing-or-single-level',
                            adjustedRand=float(adjusted_rand_score(labels,values)) if available else None)

result=dict(status='passed' if gap<=0.02 else 'reference-quality-gap',reportSHA256=hashlib.sha256(a.report.read_bytes()).hexdigest(),
 cells=n,clusters=len(partition),nativeModularity=c['modularity'],independentModularity=objective,disconnectedCommunities=disconnected,
 referenceMaximumModularityGap=gap,references=references,metadataAssociations=associations,
 qualification='Objective and heuristic reference comparison; not ground-truth cell annotation, donor integration or biological validation',
 versions={name:importlib.metadata.version(name) for name in ['numpy','scipy','networkx','scikit-learn']})
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(result,indent=2))
assert gap<=0.02, 'Native objective materially below the reference runs; retain report and investigate'
