#!/usr/bin/env python3
"""Verify native annotation preservation and full source-guide pseudobulk counts."""
import argparse
import hashlib
import json
from pathlib import Path
import anndata
import h5py
import numpy as np
from scipy import sparse
from prepare_ingestion import sha, write, SOURCE, COLUMN

def equal(left, right):
    left, right = np.asarray(left), np.asarray(right)
    assert left.shape == right.shape
    if left.dtype.kind in 'fc': assert np.array_equal(left, right, equal_nan=True)
    else: assert np.array_equal(left, right)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('source','annotated','annotation_receipt','prepared','pseudobulk','out'):
        p.add_argument('--'+key.replace('_','-'),type=Path,required=True)
    a=p.parse_args()
    assert not a.out.exists()
    assert sha(a.source)==SOURCE
    receipt=json.loads(a.annotation_receipt.read_text())
    assert bytes(receipt['source']['bytes']).hex()==SOURCE
    assert bytes(receipt['output']['bytes']).hex()==sha(a.annotated)
    reference=np.load(a.prepared/'reference.npz')
    identity=json.loads((a.prepared/'identities.json').read_text())
    audited=json.loads((a.prepared/'audit.json').read_text())
    checked=[]
    with h5py.File(a.source) as before,h5py.File(a.annotated) as after:
        def visit(path):
            b,c=before[path],after[path]
            assert isinstance(b,h5py.Dataset)==isinstance(c,h5py.Dataset)
            for name,value in b.attrs.items():
                if path=='obs' and name=='column-order':
                    equal(list(value)+[COLUMN],c.attrs[name])
                else: equal(value,c.attrs[name])
            assert set(b.attrs)==set(c.attrs)
            if isinstance(b,h5py.Group):
                expected=set(b)
                if path=='obs': expected.add(COLUMN)
                if path=='uns': expected.add('numivivo_edits')
                assert set(c)==expected
            if isinstance(b,h5py.Dataset):
                assert b.dtype==c.dtype and b.shape==c.shape
                assert (b.chunks,b.compression,b.compression_opts,b.shuffle,b.fletcher32)==(c.chunks,c.compression,c.compression_opts,c.shuffle,c.fletcher32)
                if b.ndim==0: equal(b[()],c[()])
                else:
                    for start in range(0,len(b),65536): equal(b[start:start+65536],c[start:start+65536])
                checked.append(path)
        visit('/')
        before.visit(visit)
        equal(after['obs/'+COLUMN+'/codes'][:],reference['derivedCodes'])
        equal(before['obs/perturbation/codes'][:],reference['sourceCodes'])
        journal=after['uns/numivivo_edits']
        plan_id=bytes(receipt['plan']['bytes']).hex()
        assert list(journal)==[plan_id]
        recorded=json.loads(journal[plan_id+'/plan_json'].asstr()[()])
        assert recorded==json.loads((a.prepared/'annotation-plan.json').read_text())
    obj=anndata.read_h5ad(a.annotated,backed='r')
    try:
        assert obj.shape==(65337,32738)
        assert obj.obs['perturbation'].isna().sum()==2613
        assert not obj.obs[COLUMN].isna().any()
        assert obj.obs.index.tolist()==identity['barcodes']
        assert obj.var['ensembl_id'].tolist()==identity['features']
    finally: obj.file.close()
    report=json.loads((a.pseudobulk/'report.json').read_text())
    bulk=report['pseudobulk']; m=bulk['matrix']; groups=bulk['groups']
    assert [g['condition'] for g in groups]==identity['groups']
    assert bulk['featureIDs']==identity['features']
    actual=sparse.csr_matrix((np.array(m['counts'],dtype=np.int64),m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount']))
    equal(actual.toarray(),reference['counts'])
    assert report['canonicalNonzeros']==audited['canonicalNonzeros']
    assert len(report['quality'])==65337
    for name in ['totalCounts','detectedFeatures','mitochondrialCounts']:
        equal([q[name] for q in report['quality']],reference[name])
    assert [c['barcode'] for c in report['metadata']['cells']]==identity['barcodes']
    for g in groups:
        assert len(g['sourceCellIndices'])==audited['groupCellCounts'][g['condition']]
        assert g['biologicalReplicateID']=='K562-pooled-replication-unresolved'
        assert 'donorID' not in g
        expected=np.flatnonzero(obj.obs[COLUMN].to_numpy()==g['condition']).tolist()
        assert g['sourceCellIndices']==expected
    assert sorted(i for g in groups for i in g['sourceCellIndices'])==list(range(65337))
    write(a.out,dict(status='passed',sourceSHA256=SOURCE,annotatedSHA256=sha(a.annotated),
          reportSHA256=sha(a.pseudobulk/'report.json'),originalDatasetsExactlyPreserved=checked,
          cells=65337,features=32738,sourceEntries=audited['canonicalNonzeros'],groups=len(groups),
          aggregateCountsExact=True,cellQualityCountsExact=True,sourceCellMembershipExact=True,
          annDataBackedRead=True,missingAssignmentsRetained=2613,controlsVerified=False,
          predictionQualification=False))
    print(a.out.read_text())

if __name__=='__main__': main()
