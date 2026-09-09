#!/usr/bin/env python3
"""Convert all six audited author UMI files to CSC H5AD, one gene at a time.

No dense cells-by-genes matrix, filtering, inferred donor pairing or DE claim.
Independent exact cell totals/detected genes and sample sums are retained.
"""
import argparse,gzip,hashlib,json
from contextlib import ExitStack
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-dir',type=Path,required=True)
p.add_argument('--audit',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();audit=json.loads(a.audit.read_text());records=audit['files']
assert len(records)==6 and audit['cells']==13863 and audit['nonzeros']==32848185
a.out.mkdir(parents=True,exist_ok=False)
for record in records:
    with (a.source_dir/record['file']).open('rb') as f:
        assert hashlib.file_digest(f,'sha256').hexdigest()==record['sha256']
samples=[];barcodes=[];sample_ids=[];genes=[]
with ExitStack() as stack:
    inputs=[stack.enter_context(gzip.open(a.source_dir/r['file'],'rt')) for r in records]
    for record,f in zip(records,inputs):
        sid=record['file'].split('_filtered')[0];bc=f.readline().split()
        assert len(bc)==record['cells'] and len(set(bc))==len(bc)
        samples.append(dict(id=sid,biologicalReplicateID=sid,condition='unstimulated' if '_unst' in sid else 'LPS6',batchID='unreported',organism='NCBITaxon:10090'))
        barcodes.extend(bc);sample_ids.extend([sid]*len(bc))
    offsets=np.cumsum([0]+[r['cells'] for r in records]);n=int(offsets[-1]);m=records[0]['features']
    totals=np.zeros(n,dtype=np.uint64);detected=np.zeros(n,dtype=np.int64);bulk=np.zeros((6,m),dtype=np.uint64)
    # AnnData owns the metadata encoding; the count payload is written by column.
    obs=pd.DataFrame(dict(sample=sample_ids,barcode=barcodes),index=[s+'__'+b for s,b in zip(sample_ids,barcodes)])
    path=a.out/'prepared.h5ad'
    ad.AnnData(X=sparse.csc_matrix((n,m),dtype=np.uint64),obs=obs,var=pd.DataFrame(index=[str(i) for i in range(m)])).write_h5ad(path)
    with h5py.File(path,'r+') as h:
        del h['X'];g=h.create_group('X');g.attrs['encoding-type']='csc_matrix';g.attrs['encoding-version']='0.1.0';g.attrs['shape']=[n,m]
        data=g.create_dataset('data',shape=(0,),maxshape=(None,),dtype='uint64',chunks=(65536,),compression='gzip',compression_opts=1)
        indices=g.create_dataset('indices',shape=(0,),maxshape=(None,),dtype='int32',chunks=(65536,),compression='gzip',compression_opts=1)
        pointers=[0]
        for j in range(m):
            pieces=[];rows=[];gene=None
            for k,f in enumerate(inputs):
                identifier,values=f.readline().split(maxsplit=1)
                if gene is None:gene=identifier
                assert gene==identifier
                raw=np.fromstring(values,sep=' ',dtype=np.float64)
                assert len(raw)==records[k]['cells'] and np.isfinite(raw).all() and (raw>=0).all() and (raw<=1_000_000).all() and (raw==np.floor(raw)).all()
                counts=raw.astype(np.uint64);active=np.flatnonzero(counts)
                lo,hi=offsets[k:k+2];totals[lo:hi]+=counts;detected[lo:hi]+=(counts>0)
                bulk[k,j]=counts.sum(dtype=np.uint64);pieces.append(counts[active]);rows.append((active+lo).astype(np.int32))
            genes.append(gene);values=np.concatenate(pieces);row_indices=np.concatenate(rows)
            start=pointers[-1];end=start+len(values);data.resize((end,));indices.resize((end,))
            data[start:end]=values;indices[start:end]=row_indices;pointers.append(end)
        assert all(not f.readline() for f in inputs)
        assert len(set(genes))==m and pointers[-1]==audit['nonzeros']
        assert hashlib.sha256(('\n'.join(genes)+'\n').encode()).hexdigest()==records[0]['featureAxisSHA256']
        g.create_dataset('indptr',data=np.array(pointers,dtype=np.int32))
        del h['var/_index'];h['var'].create_dataset('_index',data=np.array(genes,dtype=object),dtype=h5py.string_dtype())
        h['var/_index'].attrs['encoding-type']='string-array';h['var/_index'].attrs['encoding-version']='0.2.0'
    assert [int(x) for x in bulk.sum(axis=1)]==[r['totalUMIs'] for r in records]
    np.savez_compressed(a.out/'reference.npz',totals=totals,detected=detected,bulk=bulk,genes=np.array(genes),sample_ids=np.array([s['id'] for s in samples]))
mapping=dict(schemaVersion=1,id='hagai-original-mouse-all',evidence='measured',sourceDescription='E-MTAB-6754 author QC cluster0, six mouse source samples; donor pairing not yet qualified',countUnit='umiCount',matrixPath='X',samples=samples,sampleColumn='sample',barcodeColumn='barcode',mitochondrialFeatureIDs=[])
(a.out/'plan.json').write_text(json.dumps(dict(schemaVersion=1,mapping=mapping,contrasts=[]),indent=2)+'\n')
with path.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
(a.out/'preparation.json').write_text(json.dumps(dict(sourceAudit=audit,h5adSHA256=sha,cells=n,features=m,nonzeros=pointers[-1],scope='All six audited source files; count/QC only; no inferred donor mapping'),indent=2)+'\n')
print(json.dumps(dict(status='prepared',cells=n,features=m,nonzeros=pointers[-1],h5adSHA256=sha)),flush=True)
