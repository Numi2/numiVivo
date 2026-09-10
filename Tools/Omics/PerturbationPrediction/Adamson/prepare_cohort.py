#!/usr/bin/env python3
"""Freeze the original author's unambiguous-cell rule; no control inference or fitting."""
import argparse
from collections import Counter
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from prepare_ingestion import sha,write
from restore_geo_identities import BARCODES,GUIDES

RESTORED='a33bab097da97dc46414473d1761675940c59df8a53ee8517f60f98f950c334e'
AUTHOR='556cd498d006ac5c1553f3ef4aaa74cee4bad6444142b0c3c7553dff4ea1d7c1'
AUTHOR_URL='https://github.com/thomasmaxwellnorman/perturbseq_demo/blob/4328ac9ed29a6afe7162e1822c2f7632448667ba/perturbseq/cell_population.py#L203'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','barcodes','guides','author_code','mapping','out']:
        p.add_argument('--'+key.replace('_','-'),type=Path,required=True)
    a=p.parse_args()
    for path,expected in [(a.source,RESTORED),(a.barcodes,BARCODES),(a.guides,GUIDES),(a.author_code,AUTHOR)]:
        assert sha(path)==expected, str(path)
    a.out.mkdir(parents=True,exist_ok=False)
    barcodes=pd.read_csv(a.barcodes,header=None,sep='\t')[0]
    records=pd.read_csv(a.guides,index_col=0).reindex(barcodes).reset_index(drop=True)
    labels=records['guide identity']
    reasons=[]
    for row in records.to_dict('records'):
        failed=[]
        if pd.isna(row['guide identity']): failed.append('missing-original-guide-record')
        else:
            if row['number of cells']!=1: failed.append('not-one-assigned-identity')
            if row['good coverage']!=True: failed.append('not-good-guide-coverage')
            if row['guide identity']=='*': failed.append('ambiguous-guide-star')
        reasons.append(failed)
    indices=[i for i,failed in enumerate(reasons) if not failed]
    mask=records['number of cells'].eq(1)&records['good coverage'].eq(True)&labels.notna()&labels.ne('*')
    assert indices==np.flatnonzero(mask).tolist()
    selected_labels=labels.iloc[indices].tolist();groups=sorted(set(selected_labels))
    group_index={g:i for i,g in enumerate(groups)}
    rows=np.array([group_index[x] for x in selected_labels])
    membership=sparse.csr_matrix((np.ones(len(indices),dtype=np.int64),(rows,indices)),shape=(len(groups),len(barcodes)))
    counts=np.zeros((len(groups),32738),dtype=np.int64)
    with h5py.File(a.source) as f:
        assert f['obs/numivivo_geo_barcode'].asstr()[:].tolist()==barcodes.tolist()
        x=f['X'];assert tuple(x.attrs['shape'])==(len(barcodes),32738)
        ptr=x['indptr'][:].astype(np.int64)
        for start in range(0,32738,64):
            stop=min(start+64,32738);lo,hi=ptr[start],ptr[stop]
            data=x['data'][lo:hi]
            assert np.isfinite(data).all() and (data>=0).all() and (data==np.floor(data)).all()
            block=sparse.csc_matrix((data.astype(np.int64),x['indices'][lo:hi],ptr[start:stop+1]-lo),shape=(len(barcodes),stop-start))
            counts[:,start:stop]=(membership@block).toarray()
        features=f['var/ensembl_id'].asstr()[:].tolist()
    np.savez_compressed(a.out/'reference.npz',counts=counts,selectedRows=indices,groupForRow=rows)
    protocol=dict(rule='number of cells == 1 AND good coverage == True AND guide identity != *; missing records excluded',
        authorCodeURL=AUTHOR_URL,authorCodeSHA256=AUTHOR,controlsVerified=False,predictionsFitted=False,
        expressionFiltering=False,coverageThresholdRetuned=False)
    plan=json.loads(a.mapping.read_text())
    plan['cellSelection']=dict(source=dict(bytes=list(bytes.fromhex(RESTORED))),observationIndices=indices,
        provenance=protocol['rule']+'. Original author code '+AUTHOR_URL+' SHA256='+AUTHOR+
        '; original GEO guide CSV SHA256='+GUIDES+'. Control definitions remain unresolved; no fitted/scored predictions.')
    write(a.out/'plan.json',plan)
    write(a.out/'cohort.json',dict(sourceSHA256=RESTORED,protocol=protocol,selectedRows=indices,
        barcodes=barcodes.iloc[indices].tolist(),groups=groups,features=features,
        excluded=[dict(sourceRow=i,barcode=barcodes[i],guide=None if pd.isna(labels[i]) else labels[i],reasons=r) for i,r in enumerate(reasons) if r]))
    before=Counter(labels.dropna());after=Counter(selected_labels)
    audit=dict(sourceCells=len(barcodes),selectedCells=len(indices),excludedCells=len(barcodes)-len(indices),
        allFailureReasons=dict(Counter(v for rr in reasons for v in rr)),
        exclusiveReasonCombinations=dict(Counter('|'.join(rr) for rr in reasons if rr)),
        selectedUMIs=int(counts.sum()),selectedGuideGroups=len(groups),
        guides=[dict(guide=g,sourceCells=before[g],selectedCells=after[g]) for g in sorted(before)],
        zeroSelectedGuides=sorted(set(before)-set(after)),
        retainedGeneLikePrefixes=sorted({g.split('_')[0] for g in groups if '_pDS' in g or '_pBA' in g and not g.startswith(('62(mod)','63(mod)','Gal4-4(mod)'))}),
        protocol=protocol,sourceSHA256=RESTORED,planSHA256=sha(a.out/'plan.json'),referenceSHA256=sha(a.out/'reference.npz'))
    write(a.out/'audit.json',audit)
    print(json.dumps({k:v for k,v in audit.items() if k not in ['guides','protocol','retainedGeneLikePrefixes']},sort_keys=True))

if __name__=='__main__':main()
