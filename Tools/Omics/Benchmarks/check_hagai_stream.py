#!/usr/bin/env python3
"""Verify full-scope native streaming QC/pseudobulks against source and Scanpy."""
import argparse,hashlib,json
from importlib.metadata import version
from pathlib import Path
import anndata as ad
import numpy as np
import scanpy as sc

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--prepared',type=Path,required=True)
p.add_argument('--bundle',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
report=json.loads((a.bundle/'report.json').read_text())
receipt=json.loads((a.bundle/'receipt.json').read_text())
reference=np.load(a.prepared/'reference.npz')
obj=ad.read_h5ad(a.prepared/'prepared.h5ad')
assert obj.shape==(13863,22048) and obj.X.nnz==32848185
qc,_=sc.pp.calculate_qc_metrics(obj,percent_top=None,log1p=False,inplace=False)
native_totals=np.array([q['totalCounts'] for q in report['quality']],dtype=np.uint64)
native_detected=np.array([q['detectedFeatures'] for q in report['quality']])
assert np.array_equal(native_totals,reference['totals']) and np.array_equal(native_totals,qc['total_counts'].to_numpy())
assert np.array_equal(native_detected,reference['detected']) and np.array_equal(native_detected,qc['n_genes_by_counts'].to_numpy())
assert all(q.get('mitochondrialFraction') is None for q in report['quality'])
bulk=report['pseudobulk'];matrix=bulk['matrix'];groups=bulk['groups']
assert len(groups)==6 and bulk['featureIDs']==reference['genes'].tolist()
members=[]
for i,group in enumerate(groups):
    assert len(group['sampleIDs'])==1 and group.get('donorID') is None
    sid=group['sampleIDs'][0];k=reference['sample_ids'].tolist().index(sid)
    lo,hi=matrix['rowOffsets'][i:i+2]
    row=np.zeros(obj.n_vars,dtype=np.uint64)
    row[np.array(matrix['featureIndices'][lo:hi],dtype=np.int64)]=matrix['counts'][lo:hi]
    assert np.array_equal(row,reference['bulk'][k])
    mask=obj.obs['sample'].to_numpy()==sid
    assert np.array_equal(row,np.asarray(obj.X[mask].sum(axis=0)).ravel())
    assert group['sourceCellIndices']==np.flatnonzero(mask).tolist()
    members.extend(group['sourceCellIndices'])
assert sorted(members)==list(range(obj.n_obs))
assert report['canonicalNonzeros']==obj.X.nnz and report['contrasts']==[]
assert [c['barcode'] for c in report['metadata']['cells']]==obj.obs['barcode'].tolist()
with (a.bundle/'original.h5ad').open('rb') as f:sha=hashlib.file_digest(f,'sha256').digest()
assert list(sha)==receipt['source']['bytes']
result=dict(status='passed-full-source-count-QC-pseudobulk-comparison',cells=obj.n_obs,features=obj.n_vars,nonzeros=obj.X.nnz,
    totalUMIs=int(native_totals.sum()),pseudobulks=len(groups),h5adSHA256=sha.hex(),receipt=receipt,
    checks=['all cell totals equal original source and Scanpy','all detected-feature counts equal original source and Scanpy','all six pseudobulks equal original source and sparse AnnData sums','all cell identities and membership retained','all ordered feature identities retained','missing mitochondrial annotations remain missing'],
    versions={name:version(name) for name in ['scanpy','anndata','h5py','numpy','scipy']},
    limitations=['No verified individual donor pairing or DE biology claim','HDF5 source is chunked; metadata and pseudobulk report remain resident','No million-cell or general out-of-core qualification'])
(a.out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','cells','features','nonzeros','totalUMIs']}))
