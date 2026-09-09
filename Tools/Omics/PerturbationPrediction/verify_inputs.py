#!/usr/bin/env python3
"""Reaggregate archived experimental H5AD counts without any cell-by-gene dense array."""
import argparse,hashlib,json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
meta=json.loads((a.input/'input.json').read_text());source=Path(meta['sourceH5AD'])
with source.open('rb') as f: fingerprint=hashlib.file_digest(f,'sha256').hexdigest()
assert fingerprint==meta['sourceH5ADSHA256']
mapping=json.loads((source.parent/'plan.json').read_text())['mapping'];assert mapping['matrixPath']=='X'
obj=ad.read_h5ad(source);assert sparse.issparse(obj.X)
raw=obj.X.tocsr();assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.array_equal(raw.data,np.floor(raw.data))
counts=pd.read_csv(a.input/'counts.tsv',sep='\t',index_col=0);samples=pd.read_csv(a.input/'samples.tsv',sep='\t').set_index('sampleID')
ids=obj.var[mapping['featureIDColumn']].astype(str).tolist() if mapping.get('featureIDColumn') else obj.var_names.tolist()
assert ids==counts.index.tolist() and counts.columns.tolist()==samples.index.tolist()
lookup={s:i for i,s in enumerate(samples.index)}
rows=np.asarray([lookup[str(s)] for s in obj.obs[mapping['sampleColumn']]],dtype=np.int64)
aggregate=sparse.csr_matrix((np.ones(len(rows),dtype=np.int64),(rows,np.arange(len(rows)))),shape=(len(samples),len(rows)))
actual=(aggregate@raw).toarray();np.testing.assert_array_equal(actual,counts.to_numpy().T)
np.testing.assert_array_equal(actual.sum(axis=1),samples.libraryCounts.to_numpy())
spec={s['id']:s for s in mapping['samples']}
for sid,row in samples.iterrows():
 assert row.donor==spec[sid]['donorID'] and row.condition==spec[sid]['condition']
a.out.write_text(json.dumps({'passed':True,'sourceSHA256':fingerprint,'cells':obj.n_obs,'features':obj.n_vars,'pseudobulks':len(samples),'sourceNonzeros':int(raw.nnz),'exactCountsAxesLibrariesDonorsConditions':True},indent=2)+'\n')
