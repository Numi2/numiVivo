#!/usr/bin/env python3
"""Audit every raw entry and independently aggregate the complete source in bounded CSR blocks."""
import argparse,hashlib,json
from pathlib import Path
import h5py,numpy as np
from anndata.io import read_elem
from scipy import sparse

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
a=p.parse_args();root=a.root;receipt=json.loads((root/'source-receipt.json').read_text())
h=hashlib.sha256()
with (root/'source.h5ad').open('rb') as f:
    for b in iter(lambda:f.read(8388608),b''):h.update(b)
assert h.hexdigest()==receipt['sha256']
with h5py.File(root/'source.h5ad','r') as f:
    obs=read_elem(f['obs']);var=read_elem(f['raw/var']);x=f['raw/X']
    assert x.attrs['encoding-type']=='csr_matrix' and tuple(x.attrs['shape'])==(len(obs),len(var))
    donors=sorted(obs.donor_id.astype(str).unique());positions={d:i for i,d in enumerate(donors)}
    sample=np.array([positions[d] for d in obs.donor_id.astype(str)],dtype=np.int64)
    offsets=x['indptr'][:];assert offsets[0]==0 and np.all(np.diff(offsets)>=0)
    assert offsets[-1]==len(x['data'])==len(x['indices'])
    aggregate=np.zeros((len(donors),len(var)),dtype=np.uint64)
    totals=np.zeros(len(obs),dtype=np.uint64);detected=np.zeros(len(obs),dtype=np.int64);mito=np.zeros(len(obs),dtype=np.uint64)
    mt=np.asarray(var.mito,dtype=bool);nonzeros=0;maximum=0
    for start in range(0,len(obs),512):
        stop=min(start+512,len(obs));lo=int(offsets[start]);hi=int(offsets[stop])
        raw=x['data'][lo:hi];col=x['indices'][lo:hi];ptr=offsets[start:stop+1]-lo
        assert np.isfinite(raw).all() and (raw>=0).all() and (raw==np.floor(raw)).all() and (raw<2**53).all()
        assert (col>=0).all() and (col<len(var)).all()
        values=raw.astype(np.uint64);matrix=sparse.csr_matrix((values,col,ptr),shape=(stop-start,len(var)))
        # Source is canonical CSR; reject silently duplicated coordinates.
        assert matrix.has_canonical_format
        matrix.eliminate_zeros();nonzeros+=matrix.nnz
        maximum=max(maximum,int(values.max(initial=0)))
        totals[start:stop]=np.asarray(matrix.sum(axis=1,dtype=np.uint64)).ravel()
        detected[start:stop]=np.diff(matrix.indptr)
        mito[start:stop]=np.asarray(matrix[:,mt].sum(axis=1,dtype=np.uint64)).ravel()
        assignment=sparse.csr_matrix((np.ones(stop-start,dtype=np.uint64),(sample[start:stop],np.arange(stop-start))),shape=(len(donors),stop-start))
        aggregate+=(assignment@matrix).toarray() # donor × gene only
    assert np.array_equal(aggregate.sum(axis=1),np.bincount(sample,weights=totals,minlength=len(donors)).astype(np.uint64))
    assert sum(int(v) for v in totals)==sum(int(v) for v in aggregate.ravel())
    np.savez_compressed(root/'reference-counts.npz',counts=aggregate,donors=np.array(donors),features=np.array(var.index.astype(str)))
    np.savez_compressed(root/'reference-qc.npz',totals=totals,detected=detected,mitochondrialCounts=mito,barcodes=np.array(obs.index.astype(str)),sample=sample)
report=dict(status='all-source-counts-audited',sourceSHA256=receipt['sha256'],cells=len(obs),features=len(var),donors=len(donors),
    storedEntries=int(offsets[-1]),canonicalNonzeros=nonzeros,maximumStoredCount=maximum,totalUMIs=sum(int(v) for v in totals),
    aggregateNonzeros=int(np.count_nonzero(aggregate)),allCountsIntegral=True,cellBlockRows=512,
    referenceCountsSHA256=hashlib.sha256((root/'reference-counts.npz').read_bytes()).hexdigest(),
    referenceQCSHA256=hashlib.sha256((root/'reference-qc.npz').read_bytes()).hexdigest())
(root/'count-audit.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(json.dumps(report,indent=2))
