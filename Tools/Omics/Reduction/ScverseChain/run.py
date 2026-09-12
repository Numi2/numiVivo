from pathlib import Path
import json,time,hashlib
import numpy as np,anndata as ad,scanpy as sc
from scipy import sparse
r=Path('/Users/n/numivivo-scverse-chain-20260912');source=Path('/Users/n/numivivo-pca-borrowed-20260912/pca-1/original.h5ad')
protocol={'cohort':'Kang','cells':24673,'normalizationTarget':10000,'HVG':{'flavor':'seurat','n_top_genes':2000,'n_bins':20},'PCA':{'n_comps':20,'svd_solver':'arpack','dtype':'float64','random_state':7},'neighbors':{'n_neighbors':20,'n_pcs':20,'use_rep':'X_pca','knn':True,'method':'umap','metric':'euclidean','transformer':'sklearn','random_state':7},'gate':'exact nonself neighbor sets and graph edge coordinates; weight and distance differences retained','scope':'independent raw-count-to-graph numerical comparison, single run, not biological validation or matched timing'}
(r/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
t=time.monotonic();a=ad.read_h5ad(source);assert sparse.issparse(a.X);a.X=a.X.astype(np.float64);sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,**protocol['HVG']);sc.pp.pca(a,**protocol['PCA']);sc.settings.n_jobs=1;sc.pp.neighbors(a,**protocol['neighbors']);elapsed=time.monotonic()-t
for name in ['distances','connectivities']:sparse.save_npz(r/(name+'.npz'),a.obsp[name])
np.save(r/'scores.npy',a.obsm['X_pca']);(r/'selected-features.json').write_text(json.dumps(a.var_names[a.var['highly_variable']].tolist()))
dist=a.obsp['distances'].tocsr();conn=a.obsp['connectivities'].tocsr();conn.sort_indices();results={}
for mode,folder in [('cpu','0-cpu'),('metal','2-metal')]:
 path=Path('/Users/n/numivivo-metal-scverse-20260912')/folder/'graph.json';g=json.loads(path.read_text());n=len(g['cells']);assert n==a.n_obs
 ids=np.asarray(g['neighborIndices']).reshape(n,20);values=np.asarray(g['neighborDistances']).reshape(n,20);missing=extra=0;error=0.;entries=0
 for i in range(n):
  lo,hi=dist.indptr[i:i+2];ref={int(j):float(v) for j,v in zip(dist.indices[lo:hi],dist.data[lo:hi]) if j!=i};candidate=dict(zip(ids[i,1:].tolist(),values[i,1:].tolist()));missing+=len(ref.keys()-candidate.keys());extra+=len(candidate.keys()-ref.keys());entries+=len(ref)
  if ref.keys()&candidate.keys():error=max(error,max(abs(ref[j]-candidate[j]) for j in ref.keys()&candidate.keys()))
 graph=sparse.csr_matrix((g['weights'],g['columnIndices'],g['rowOffsets']),shape=conn.shape);graph.sort_indices();same=np.array_equal(graph.indptr,conn.indptr) and np.array_equal(graph.indices,conn.indices)
 results[mode]={'nativeGraphSHA256':hashlib.sha256(path.read_bytes()).hexdigest(),'referenceNeighborEntries':entries,'missingNeighbors':missing,'extraNeighbors':extra,'maximumCommonNeighborDistanceDifference':error,'nativeEdges':int(graph.nnz),'referenceEdges':int(conn.nnz),'exactEdgeCoordinates':same,'maximumWeightDifference':float(np.max(np.abs(graph.data-conn.data))) if same else None,'gatePass':missing==extra==0 and same}
(r/'results.json').write_text(json.dumps({'comparisons':results,'scanpyReadNormalizeHVGPCAGraphSeconds':elapsed,'sourceSHA256':hashlib.sha256(source.read_bytes()).hexdigest()},indent=2)+'\n');print(json.dumps(results))
