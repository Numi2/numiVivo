#!/usr/bin/env python3
"""Restore GEO barcode/guide identity without replacing deposited AnnData fields."""
import argparse
from collections import Counter
import json
from pathlib import Path
import anndata
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from prepare_ingestion import SOURCE,sha,write

BARCODES='1e0d820343d0c6e17bdab4fef96f4446e223b94861942ccf3d7225818e009836'
GUIDES='8b40be7a2280c1713bf5a1eb828aad46e70ad3a7000b5f5a1322e51e06c4cf7f'
MISSING='numivivo-unassigned-geo-guide'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('source','barcodes','guides','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();assert sha(a.source)==SOURCE and sha(a.barcodes)==BARCODES and sha(a.guides)==GUIDES
    a.out.mkdir(parents=True,exist_ok=False)
    barcodes=pd.read_csv(a.barcodes,header=None,sep='\t')[0]
    guides=pd.read_csv(a.guides,index_col=0)
    assert len(barcodes)==65337 and barcodes.is_unique and guides.index.is_unique
    assert set(guides.index).issubset(set(barcodes))
    original=guides.reindex(barcodes).reset_index(drop=True)
    true_labels=original['guide identity']; labels=true_labels.fillna(MISSING)
    gem=barcodes.str.rsplit('-',n=1).str[1]
    assert set(gem)=={str(i) for i in range(1,11)}
    categories=sorted(labels.unique()); label_index={v:i for i,v in enumerate(categories)}
    rows=np.array([label_index[x] for x in labels])
    samples=labels+'@GEM-'+gem; sample_categories=sorted(samples.unique())
    sample_index={v:i for i,v in enumerate(sample_categories)}
    sample_records=[]
    for sample in sample_categories:
        i=int(np.flatnonzero(samples.to_numpy()==sample)[0])
        sample_records.append(dict(id=sample,condition=labels[i],batchID='GSM2406681-GEM-'+gem[i],
            biologicalReplicateID='K562-pooled-replication-unresolved',organism='NCBITaxon:9606'))
    with h5py.File(a.source) as f:
        source_barcodes=f['obs/cell_barcode'].asstr()[:].tolist()
        # The deposited preprocessing strips suffixes and then makes names unique.
        reconstructed=anndata.utils.make_index_unique(pd.Index(barcodes.str.split('-').str[0]))
        assert reconstructed.tolist()==source_barcodes
        old_categories=f['obs/perturbation/categories'].asstr()[:];old_codes=f['obs/perturbation/codes'][:]
        old=pd.Series([old_categories[i] if i>=0 else None for i in old_codes])
        # Independently reproduce the curator's actual first-duplicate join.
        bare=guides.copy();bare.index=bare.index.str.split('-').str[0];bare=bare[~bare.index.duplicated()]
        replay=bare.reindex(source_barcodes)['guide identity'].reset_index(drop=True)
        assert replay.fillna('<missing>').tolist()==old.fillna('<missing>').tolist()
        changed=old.fillna('<missing>')!=true_labels.fillna('<missing>')
        restored=int((old.isna()&true_labels.notna()).sum());removed=int((old.notna()&true_labels.isna()).sum())
        different=int((changed&old.notna()&true_labels.notna()).sum())
        ids=f['var/ensembl_id'].asstr()[:].tolist();names=f['var/gene_symbol'].asstr()[:].tolist()
        mito=[ids[i] for i,n in enumerate(names) if n.startswith('MT-')]
        x=f['X'];assert tuple(x.attrs['shape'])==(65337,32738)
        membership=sparse.csr_matrix((np.ones(len(rows),dtype=np.int64),(rows,np.arange(len(rows)))),shape=(len(categories),len(rows)))
        counts=np.zeros((len(categories),32738),dtype=np.int64);pointers=x['indptr'][:].astype(np.int64)
        for start in range(0,32738,64):
            stop=min(32738,start+64);lo,hi=pointers[start],pointers[stop]
            data=x['data'][lo:hi];assert np.isfinite(data).all() and np.all(data>=0) and np.all(data==np.floor(data))
            block=sparse.csc_matrix((data.astype(np.int64),x['indices'][lo:hi],pointers[start:stop+1]-lo),shape=(65337,stop-start))
            counts[:,start:stop]=(membership@block).toarray()
        np.savez_compressed(a.out/'reference.npz',counts=counts,groupForRow=rows)
    write(a.out/'identities.json',dict(groups=categories,features=ids,names=names,barcodes=barcodes.tolist(),
        depositedBarcodes=source_barcodes,sourceGuides=[None if pd.isna(v) else v for v in true_labels],
        gemGroups=gem.tolist(),samples=samples.tolist()))
    edits=[]
    def edit(name,value):edits.append(dict(path='obs/numivivo_geo_'+name,mode='add',value=value))
    edit('barcode',dict(string=dict(shape=[len(barcodes)],values=barcodes.tolist())))
    edit('guide',dict(categorical=dict(codes=rows.tolist(),categories=categories,ordered=False)))
    edit('sample',dict(categorical=dict(codes=[sample_index[x] for x in samples],categories=sample_categories,ordered=False)))
    edit('gem_group',dict(categorical=dict(codes=[int(x)-1 for x in gem],categories=[str(i) for i in range(1,11)],ordered=False)))
    for source,name,kind,convert in [('good coverage','good_coverage','nullableBoolean',bool),
        ('number of cells','number_of_cells','nullableInt64',int),('UMI count','guide_umi_count','nullableInt64',int),
        ('read count','guide_read_count','nullableInt64',int),('coverage','guide_coverage','nullableFloat64',float)]:
        edit(name,{kind:dict(values=[None if pd.isna(x) else convert(x) for x in original[source]])})
    write(a.out/'annotation-plan.json',dict(schemaVersion=1,source=dict(bytes=list(bytes.fromhex(SOURCE))),
        provenance='Restore Adamson GSM2406681 guide records using original full GEO cell barcodes. Original barcode SHA256='+BARCODES+
        '; guide CSV SHA256='+GUIDES+'. Row order verified by exactly reconstructing the deposited barcode-renaming and erroneous metadata join. All original AnnData columns and counts retained; no control assignments or coverage filtering.',edits=edits))
    mapping=dict(schemaVersion=1,id='adamson2016-upr-geo-identities',evidence='measured',countUnit='umiCount',matrixPath='X',
        sourceDescription='Complete Adamson UPR source counts with full GEO barcode and guide identity restored in new columns. Ten GEM groups are technical batches, not independent biological replicates. No guide/control inference or coverage filtering.',
        sampleColumn='numivivo_geo_sample',barcodeColumn='numivivo_geo_barcode',featureIDColumn='ensembl_id',groupColumn='cell_line',
        samples=sample_records,mitochondrialFeatureIDs=mito)
    write(a.out/'pseudobulk-plan.json',dict(schemaVersion=1,mapping=mapping,contrasts=[]))
    write(a.out/'identity-audit.json',dict(sourceSHA256=SOURCE,geoBarcodesSHA256=BARCODES,geoGuidesSHA256=GUIDES,
        cells=65337,features=32738,originalGuideRows=len(guides),depositedBarcodeRenameReproducedExactly=True,
        depositedGuideJoinReproducedExactly=True,duplicateBareBarcodes=len(barcodes)-barcodes.str.split('-').str[0].nunique(),
        changedAssignments=int(changed.sum()),restoredAssignments=restored,removedUnsupportedAssignments=removed,differentGuideAssignments=different,
        missingDeposited=int(old.isna().sum()),missingOriginal=int(true_labels.isna().sum()),gemGroups=dict(Counter(gem)),
        groupCellCounts=dict(Counter(labels)),sampleCount=len(sample_records),totalCounts=int(counts.sum()),
        sourceGoodCoverageCounts={str(k):int(v) for k,v in original['good coverage'].value_counts(dropna=False).items()},
        controlsVerified=False,coverageFilteringPerformed=False,predictionFittingPerformed=False))
    print((a.out/'identity-audit.json').read_text())

if __name__=='__main__':main()
