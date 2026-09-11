#!/usr/bin/env python3
"""Full AnnData/SciPy count checks and source-bound confident-cohort preparation."""
import argparse,hashlib,importlib.metadata,json
from pathlib import Path
import anndata as ad,numpy as np
from scipy import sparse
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();r=a.root;prepared=r/'prepared-complete';reference=np.load(prepared/'reference.npz');inv=json.loads((prepared/'inventory.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
records=[]
for label,nf in [('rna',33694),('guide',64)]:
 obj=ad.read_h5ad(r/(label+'-projection/projected.h5ad'),backed='r');assert obj.shape==(40997,nf)
 try:
  records.append(dict(assay=label,shape=list(obj.shape),backedType=type(obj.X).__name__,countType=str(obj.X.dtype),observationIndexUnique=obj.obs_names.is_unique,featureIndexUnique=obj.var_names.is_unique))
  assert obj.obs_names.is_unique and obj.var_names.is_unique
  if label=='rna':
   selected=np.zeros(40997,dtype=bool);selected[reference['selected']]=True
   all_counts=np.zeros((len(inv['groups']),nf),dtype=np.int64);selected_counts=np.zeros_like(all_counts);totals=[];detected=[]
   for start in range(0,40997,512):
    stop=min(start+512,40997);x=obj.X[start:stop].astype(np.int64);x.sum_duplicates();x.eliminate_zeros()
    totals.extend(np.asarray(x.sum(axis=1)).ravel().tolist());detected.extend(np.diff(x.indptr).tolist())
    assignment=reference['group'][start:stop];cols=np.arange(stop-start)
    member=sparse.csr_matrix((np.ones(stop-start,dtype=np.int64),(assignment,cols)),shape=(len(inv['groups']),stop-start))
    all_counts+=(member@x).toarray()
    selected_member=sparse.csr_matrix((selected[start:stop].astype(np.int64),(assignment,cols)),shape=member.shape);selected_member.eliminate_zeros()
    selected_counts+=(selected_member@x).toarray()
   np.testing.assert_array_equal(all_counts,reference['counts'][:,reference['rnaMask']]);np.testing.assert_array_equal(totals,reference['rnaTotals']);np.testing.assert_array_equal(detected,reference['rnaDetected'])
   np.savez_compressed(r/'selected-reference.npz',counts=selected_counts,group=reference['group'],selected=reference['selected'],totals=np.asarray(totals),detected=np.asarray(detected),mitochondrial=reference['mitochondrial'])
 finally:obj.file.close()
plan=json.loads((prepared/'plan.json').read_text());plan['cellSelection']=dict(source=dict(bytes=list(bytes.fromhex(sha(r/'rna-projection/projected.h5ad')))),observationIndices=reference['selected'].tolist(),provenance='Original GEO full-barcode join; good_coverage=True and exact integer number_of_cells=1; controls from Replogle 2020 Pilot UPR methods; all five gemgroups separately retained; no expression-based selection')
(r/'selected-plan.json').write_text(json.dumps(plan,sort_keys=True,separators=(',',':'))+'\n')
(r/'scverse-checks.json').write_text(json.dumps(dict(status='passed',records=records,allRNAGroupCountsExact=True,allRNAQualityExact=True,selectedCells=len(reference['selected']),selectedGroups=int(np.count_nonzero(selected_counts.sum(axis=1))),versions={n:importlib.metadata.version(n) for n in ['anndata','scipy','numpy','h5py','pandas']},predictionFitted=False),indent=2)+'\n');print('Full RNA AnnData/SciPy counts exact; confident cohort prepared')
