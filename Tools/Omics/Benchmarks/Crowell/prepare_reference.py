#!/usr/bin/env python3
"""Check full native count/QC identity and prepare exact per-animal R comparisons."""
import argparse,gzip,hashlib,json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();protocol=json.loads((root/'protocol.json').read_text());assert sha(root/'counts.h5ad')==protocol['sourceH5ADSHA256']
obj=ad.read_h5ad(root/'counts.h5ad');qc=pd.read_csv(root/'scanpy-qc.tsv.gz',sep='\t',index_col=0,float_precision='round_trip');assert list(qc.index)==obj.obs_names.tolist()
records=[];feature_ids=obj.var.ENSEMBL.astype(str).tolist();symbols=obj.var.SYMBOL.astype(str).tolist();sample_conditions=dict(zip(obj.obs.sample_id.astype(str),obj.obs.group_id.astype(str)))
for case in protocol['cases']:
 folder=root/case['id'];path=folder/'native/report.json.gz'
 if not path.exists():
  records.append(dict(case=case['id'],cellGroup=case['cellGroup'],status='native-unavailable',expectedReplicationGate=case['expectedReplicationGate']));continue
 report=json.loads(gzip.decompress(path.read_bytes()));c=report['contrasts'][0];b=report['pseudobulk'];design=c['design'];assert c['request']['design']=='independentReplicates';assert c['request']['cellGroup']==case['cellGroup']
 assert report['canonicalNonzeros']==protocol['nonzeros'];assert b['featureIDs']==feature_ids
 assert [(v['id'],v['name']) for v in report['metadata']['features']]==list(zip(feature_ids,symbols))
 assert [(v['sampleID'],v['barcode'],v['group']) for v in report['metadata']['cells']]==list(zip(obj.obs.sample_id.astype(str),obj.obs.barcode.astype(str),obj.obs.cluster_id.astype(str)))
 observed=report['quality'];np.testing.assert_array_equal([v['totalCounts'] for v in observed],qc.total_counts);np.testing.assert_array_equal([v['detectedFeatures'] for v in observed],qc.n_genes_by_counts)
 np.testing.assert_array_equal([v['mitochondrialCounts'] for v in observed],qc.total_counts_mitochondrial)
 mito_count=int(obj.var.SYMBOL.str.lower().str.startswith('mt-').sum())
 assert all(v['mitochondrialFeatureCount']==mito_count for v in observed)
 if mito_count:np.testing.assert_allclose([v['mitochondrialFraction'] for v in observed],qc.pct_counts_mitochondrial/100,rtol=0,atol=2e-15)
 else:assert all('mitochondrialFraction' not in v for v in observed)
 rows=[];cols=[]
 for i,g in enumerate(b['groups']):
  assert len(g['sampleIDs'])==1 and g['donorID']==g['biologicalReplicateID']==g['sampleIDs'][0]
  sid=g['sampleIDs'][0];assert g['condition']==sample_conditions[sid]
  expected=np.flatnonzero((obj.obs.sample_id.astype(str)==sid)&(obj.obs.cluster_id.astype(str)==g['cellGroup']))
  np.testing.assert_array_equal(expected,g['sourceCellIndices']);rows.extend([i]*len(expected));cols.extend(expected.tolist())
 assert len(cols)==len(set(cols))==obj.n_obs
 aggregate=sparse.csr_matrix((np.ones(len(rows),dtype=np.int64),(rows,cols)),shape=(len(b['groups']),obj.n_obs));all_counts=(aggregate@obj.X).toarray();m=b['matrix'];native=sparse.csr_matrix((m['counts'],m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount'])).toarray();np.testing.assert_array_equal(all_counts,native)
 y=all_counts[design['sourcePseudobulkIndices']];ids=[o['sampleIDs'][0] for o in design['observations']];assert len(ids)==8 and len(set(ids))==8 and design['rows']==[[1,0 if sample_conditions[s]=='Vehicle' else 1] for s in ids] and design['contrast']==[0,1]
 np.testing.assert_array_equal(y.sum(axis=1),design['libraryCounts']);eligible=(y.sum(axis=0)>=10)&((y>0).sum(axis=0)>=3);assert np.array_equal(~eligible,[f['status']=='filteredLowExpression' for f in c['features']])
 # Independent standard median-ratio normalization uses only all-positive genes.
 reference=(y>0).all(axis=0);geomean=np.exp(np.log(y[:,reference]).mean(axis=0));sf=np.median(y[:,reference]/geomean,axis=1);sf/=np.exp(np.log(sf).mean());np.testing.assert_allclose(sf,design['sizeFactorValues'],rtol=2e-13,atol=1e-14)
 out=folder/'reference-input';out.mkdir(exist_ok=False)
 pd.DataFrame(y.T,index=feature_ids,columns=ids).to_csv(out/'counts.tsv',sep='\t',index_label='featureID')
 pd.DataFrame(design['rows'],index=ids,columns=design['columnNames']).to_csv(out/'design.tsv',sep='\t',index_label='sampleID')
 pd.DataFrame(dict(sampleID=ids,donor=ids,condition=[sample_conditions[s] for s in ids],sizeFactor=design['sizeFactorValues'],libraryCounts=design['libraryCounts'])).to_csv(out/'samples.tsv',sep='\t',index=False)
 pd.DataFrame(c['features']).assign(eligible=eligible).to_csv(out/'native.tsv',sep='\t',index=False)
 meta=dict(schemaVersion=1,contrast=design['contrast'],minimumFeatureCounts=10,minimumExpressingPseudobulks=3,expectedGenes=[gene for gene,symbol in zip(feature_ids,symbols) if symbol in protocol['expectedBiology']['predeclaredGeneSymbols']],observations=8,features=len(feature_ids),eligibleFeatures=int(eligible.sum()),nativeTrend=c['negativeBinomial']['trend']['method'],sourceReportSHA256=hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest(),sourceH5ADSHA256=sha(root/'counts.h5ad'),exactCountsAndAxes=True,design='Eight independent animals: intercept plus LPS-minus-Vehicle. No donor pairing or batch covariate.',qualification='Third public treatment study; descriptive reference comparison, not FDR or causal qualification')
 (out/'input.json').write_text(json.dumps(meta,indent=2)+'\n');records.append(dict(case=case['id'],cellGroup=case['cellGroup'],status='all-source-QC-counts-and-axes-exact',sourceCells=obj.n_obs,sourceGenes=obj.n_vars,pseudobulkGroups=len(b['groups']),caseCells=case['cells'],eligibleFeatures=int(eligible.sum()),nativeTested=c['testedFeatures'],maximumSizeFactorDifference=float(np.max(np.abs(sf-np.array(design['sizeFactorValues']))))))
(root/'preparation-check.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records,indent=2))
