#!/usr/bin/env python3
"""Check restored identities against original GEO records and native counts."""
import argparse
import json
from pathlib import Path
import anndata
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from prepare_ingestion import sha,write,SOURCE
from restore_geo_identities import BARCODES,GUIDES,MISSING

def equal(a,b):
    a,b=np.asarray(a),np.asarray(b);assert a.shape==b.shape
    assert np.array_equal(a,b,equal_nan=True) if a.dtype.kind in 'fc' else np.array_equal(a,b)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('source','barcodes','guides','restored','old_prepared','out'):
        p.add_argument('--'+key.replace('_','-'),type=Path,required=True)
    a=p.parse_args();assert not a.out.exists()
    assert sha(a.source)==SOURCE and sha(a.barcodes)==BARCODES and sha(a.guides)==GUIDES
    barcodes=pd.read_csv(a.barcodes,header=None,sep='\t')[0]
    metadata=pd.read_csv(a.guides,index_col=0).reindex(barcodes).reset_index(drop=True)
    plan=json.loads((a.restored/'annotation-plan.json').read_text())
    receipt=json.loads((a.restored/'annotation-receipt.json').read_text())
    output=a.restored/'annotated.h5ad'
    assert bytes(receipt['source']['bytes']).hex()==SOURCE and bytes(receipt['output']['bytes']).hex()==sha(output)
    checked=[]
    with h5py.File(a.source) as before,h5py.File(output) as after:
        def visit(path):
            b,c=before[path],after[path];assert isinstance(b,h5py.Dataset)==isinstance(c,h5py.Dataset)
            assert set(b.attrs)==set(c.attrs)
            for key,value in b.attrs.items():
                if path=='obs' and key=='column-order':
                    equal(list(value)+[e['path'].split('/')[1] for e in plan['edits']],c.attrs[key])
                else:equal(value,c.attrs[key])
            if isinstance(b,h5py.Group):
                extra={e['path'].split('/')[1] for e in plan['edits']} if path=='obs' else ({'numivivo_edits'} if path=='uns' else set())
                assert set(c)==set(b)|extra
            else:
                assert b.dtype==c.dtype and b.shape==c.shape
                assert (b.chunks,b.compression,b.compression_opts,b.shuffle,b.fletcher32)==(c.chunks,c.compression,c.compression_opts,c.shuffle,c.fletcher32)
                if b.ndim==0:equal(b[()],c[()])
                else:
                    for start in range(0,len(b),65536):equal(b[start:start+65536],c[start:start+65536])
                checked.append(path)
        visit('/');before.visit(visit)
        journal=after['uns/numivivo_edits'];key=bytes(receipt['plan']['bytes']).hex()
        assert list(journal)==[key]
        assert json.loads(journal[key+'/plan_json'].asstr()[()])==plan
    adata=anndata.read_h5ad(output,backed='r')
    try:
        obs=adata.obs
        assert obs['numivivo_geo_barcode'].tolist()==barcodes.tolist()
        labels=metadata['guide identity'].fillna(MISSING)
        assert obs['numivivo_geo_guide'].tolist()==labels.tolist()
        gem=barcodes.str.rsplit('-',n=1).str[-1]
        assert obs['numivivo_geo_gem_group'].tolist()==gem.tolist()
        assert obs['numivivo_geo_sample'].tolist()==(labels+'@GEM-'+gem).tolist()
        for original,derived in [('good coverage','good_coverage'),('number of cells','number_of_cells'),
            ('UMI count','guide_umi_count'),('read count','guide_read_count'),('coverage','guide_coverage')]:
            actual=obs['numivivo_geo_'+derived].reset_index(drop=True);expected=metadata[original]
            equal(actual.isna(),expected.isna());mask=~expected.isna();equal(actual[mask],expected[mask])
        new_quality=np.load(a.old_prepared/'reference.npz')
        report=json.loads((a.restored/'pseudobulk/report.json').read_text())
        native=report['pseudobulk'];matrix=native['matrix'];groups=native['groups']
        identity=json.loads((a.restored/'identities.json').read_text());reference=np.load(a.restored/'reference.npz')
        assert [g['condition'] for g in groups]==identity['groups']
        assert native['featureIDs']==adata.var['ensembl_id'].tolist()==identity['features']
        actual=sparse.csr_matrix((np.array(matrix['counts'],dtype=np.int64),matrix['featureIndices'],matrix['rowOffsets']),shape=(matrix['cellCount'],matrix['featureCount']))
        equal(actual.toarray(),reference['counts'])
        equal(actual.sum(axis=0),new_quality['counts'].sum(axis=0)[None,:])
        for column in ['totalCounts','detectedFeatures','mitochondrialCounts']:
            equal([q[column] for q in report['quality']],new_quality[column])
        assert [c['barcode'] for c in report['metadata']['cells']]==barcodes.tolist()
        assert [c['sampleID'] for c in report['metadata']['cells']]==identity['samples']
        assert len(report['metadata']['samples'])==1106
        for group in groups:
            indices=np.flatnonzero(labels.to_numpy()==group['condition']).tolist()
            assert group['sourceCellIndices']==indices
            assert group['batchIDs']==sorted({'GSM2406681-GEM-'+gem[i] for i in indices})
            assert group['sampleIDs']==sorted({identity['samples'][i] for i in indices})
            assert group['biologicalReplicateID']=='K562-pooled-replication-unresolved' and 'donorID' not in group
        assert sorted(i for g in groups for i in g['sourceCellIndices'])==list(range(len(barcodes)))
        original_labels=obs['perturbation'].reset_index(drop=True).astype(object)
        changed=original_labels.fillna(MISSING)!=labels
        reassigned_counts=int(new_quality['totalCounts'][changed.to_numpy()].sum())
    finally:adata.file.close()
    write(a.out,dict(status='passed-native-GEO-identity-restoration',sourceSHA256=SOURCE,
        annotatedSHA256=sha(output),reportSHA256=sha(a.restored/'pseudobulk/report.json'),
        allOriginalDatasetsPreserved=checked,restoredMetadataColumnsExact=9,
        cells=65337,features=32738,sourceEntries=237812947,sourceGroups=115,technicalGEMGroups=10,sampleMappings=1106,
        changedAssignments=int(changed.sum()),countsInReassignedCells=reassigned_counts,
        nativeAggregateCountsExact=True,nativePerCellQualityExact=True,perGeneGlobalTotalsUnchanged=True,
        nativeOriginalBarcodeAndBatchMembershipExact=True,controlsVerified=False,predictionQualification=False))
    print(a.out.read_text())

if __name__=='__main__':main()
