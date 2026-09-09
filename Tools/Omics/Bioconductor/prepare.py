#!/usr/bin/env python3
"""Bind independent pseudobulk counts to the exact native observations/design."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
p=argparse.ArgumentParser();p.add_argument('--report',type=Path,required=True);p.add_argument('--counts',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--expected',nargs='+',required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
r=json.loads(a.report.read_text());assert len(r['contrasts'])==1
c=r['contrasts'][0];d=c['design'];b=r['pseudobulk'];m=b['matrix']
assert c['request']['model']=='negativeBinomial'
ids=[o['sampleIDs'][0] for o in d['observations']];assert all(len(o['sampleIDs'])==1 for o in d['observations'])
reference=pd.read_csv(a.counts,sep='\t',index_col=0)
assert list(reference.index)==ids
# Independently aggregate all archived H5AD genes before applying any filter.
source=a.report.parent/'original.h5ad'
obj=ad.read_h5ad(source);mapping=json.loads((a.report.parent/'plan.json').read_text())['mapping']
assert mapping['matrixPath']=='X'
feature_ids=obj.var[mapping['featureIDColumn']].astype(str).tolist() if mapping.get('featureIDColumn') else obj.var_names.astype(str).tolist()
assert feature_ids==b['featureIDs']
barcodes=obj.obs[mapping['barcodeColumn']].astype(str).tolist() if mapping.get('barcodeColumn') else obj.obs_names.astype(str).tolist()
assert list(zip(obj.obs[mapping['sampleColumn']].astype(str),barcodes))==[(cell['sampleID'],cell['barcode']) for cell in r['metadata']['cells']]
rows=[];cols=[]
for i,o in enumerate(d['observations']):
 rows.extend([i]*len(o['sourceCellIndices']));cols.extend(o['sourceCellIndices'])
assert len(cols)==len(set(cols))
aggregate=sparse.csr_matrix((np.ones(len(rows),dtype=np.int64),(rows,cols)),shape=(len(ids),obj.n_obs))
assert len(ids)*obj.n_vars<=10_000_000
raw=obj.X.tocsr();assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.array_equal(raw.data,np.floor(raw.data))
y=(aggregate@raw).toarray();assert y.max()<=np.iinfo(np.int32).max
# Only pseudobulk is dense. Exact agreement includes the filtered-out genes.
native=np.zeros((len(ids),m['featureCount']),dtype=np.int64)
for i,source_row in enumerate(d['sourcePseudobulkIndices']):
 assert b['groups'][source_row]==d['observations'][i]
 lo,hi=m['rowOffsets'][source_row:source_row+2];native[i,np.array(m['featureIndices'][lo:hi],dtype=int)]=m['counts'][lo:hi]
np.testing.assert_array_equal(y,native);np.testing.assert_array_equal(y.sum(axis=1),d['libraryCounts'])
eligible=(y.sum(axis=0)>=c['request']['minimumFeatureCounts'])&((y>0).sum(axis=0)>=c['request']['minimumExpressingPseudobulks'])
assert np.array_equal(~eligible,[f['status']=='filteredLowExpression' for f in c['features']])
assert list(reference.columns)==np.array(feature_ids)[eligible].tolist()
np.testing.assert_array_equal(reference.to_numpy(),y[:,eligible])
counts=pd.DataFrame(y,index=ids,columns=feature_ids)
counts.T.to_csv(a.out/'counts.tsv',sep='\t',index_label='featureID')
pd.DataFrame(d['rows'],index=ids,columns=d['columnNames']).to_csv(a.out/'design.tsv',sep='\t',index_label='sampleID')
pd.DataFrame(dict(sampleID=ids,donor=[o['donorID'] for o in d['observations']],condition=[o['condition'] for o in d['observations']],sizeFactor=d['sizeFactorValues'],libraryCounts=d['libraryCounts'])).to_csv(a.out/'samples.tsv',sep='\t',index=False)
pd.DataFrame(c['features']).assign(eligible=eligible).to_csv(a.out/'native.tsv',sep='\t',index=False)
meta=dict(schemaVersion=1,contrast=d['contrast'],minimumFeatureCounts=c['request']['minimumFeatureCounts'],minimumExpressingPseudobulks=c['request']['minimumExpressingPseudobulks'],expectedGenes=a.expected,observations=len(ids),features=len(eligible),eligibleFeatures=int(eligible.sum()),nativeTrend=c['negativeBinomial']['trend']['method'],sourceReport=str(a.report.resolve()),sourceReportSHA256=hashlib.sha256(a.report.read_bytes()).hexdigest(),independentCounts=str(a.counts.resolve()),independentCountsSHA256=hashlib.sha256(a.counts.read_bytes()).hexdigest(),sourceH5AD=str(source),sourceH5ADSHA256=hashlib.sha256(source.read_bytes()).hexdigest(),exactCountsAndAxes=True,design='Exact native numeric matrix with paired donor fixed effects; no invented batch covariate',qualification='Two public study contrasts; comparison does not calibrate FDR or establish biological truth')
(a.out/'input.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
