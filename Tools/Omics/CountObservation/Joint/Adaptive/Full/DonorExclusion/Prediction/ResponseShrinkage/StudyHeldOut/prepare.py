from pathlib import Path
import json,hashlib
import numpy as np,anndata as ad
from scipy import sparse
s=Path(__file__).parent;gse=Path('/Users/home/numivivo-gse181897-20260911');src=gse/'prediction-inputs';fr=json.loads((src/'input-freeze.json').read_text());bindings={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bound(p,expected):
 assert sha(p)==expected,str(p);bindings[str(p)]=expected;return p
for name in ['panel.json','cohort.json','primary-condition-roles.json']:bound(src/name,fr['files'][name])
roles=json.loads((src/'primary-condition-roles.json').read_text());assert roles['conditionMeaningQualified'] and roles['mapping']=={'B':'IFNB','C':'control'}
panel=json.loads((src/'panel.json').read_text());assert len(panel)==11800 and len(set(panel))==len(panel);cohorts={}
for origin in ['Kang','HIRISA']:
 for name in ['training.h5ad','training.json']:bound(src/origin/name,fr['files'][origin+'/'+name])
 a=ad.read_h5ad(src/origin/'training.h5ad');plan=json.loads((src/origin/'training.json').read_text());samples={v['id']:v for v in plan['mapping']['samples']};lookup={(samples[str(v)]['donorID'],samples[str(v)]['condition']):i for i,v in enumerate(a.obs['sample'])};donors=sorted(set(k[0] for k in lookup));assert len(lookup)==2*len(donors)
 counts=a.X.toarray().astype(np.uint64);totals=counts.sum(axis=1,dtype=np.uint64);assert np.all(totals>0);idx=a.var_names.get_indexer(panel);assert np.all(idx>=0) and a.var_names.is_unique;logs=np.log1p(counts.astype('f8')/totals[:,None]*1e6)[:,idx]
 cohorts[origin]={'donors':donors,'controls':logs[[lookup[(d,'control')] for d in donors]],'treated':logs[[lookup[(d,'IFNB')] for d in donors]]}
handoff=gse/'handoff';bound(handoff/'input-freeze.json',fr['originalHandoffFreezeSHA256']);hf=json.loads((handoff/'input-freeze.json').read_text())
for name in ['reference-groups.json','reference-pseudobulk.npz']:bound(handoff/name,hf['files'][name])
execution=gse/'native-handoff-qualified/execution.json';bound(execution,fr['nativeAggregationExecutionSHA256']);report=gse/'native-handoff-qualified/aggregate/report.json';bound(report,json.loads(execution.read_text())['reportSHA256']);bulk=json.loads(report.read_text())['pseudobulk'];counts=sparse.load_npz(handoff/'reference-pseudobulk.npz').toarray().astype(np.uint64);groups=json.loads((handoff/'reference-groups.json').read_text());lookup={(v['expID'],v['conditionCode']):i for i,v in enumerate(groups)}
bm=bulk['matrix'];native=sparse.csr_matrix((bm['counts'],bm['featureIndices'],bm['rowOffsets']),shape=(bm['cellCount'],bm['featureCount'])).toarray().astype(np.uint64)
for row,g in enumerate(bulk['groups']):
 key=(g['donorID'].removeprefix('GSE181897:exp_id:'),g['condition'].removeprefix('source-code:'));np.testing.assert_array_equal(native[row],counts[lookup[key]])
ids=bulk['featureIDs'];idx=[ids.index(v.removeprefix('symbol|')) for v in panel];assert len(set(idx))==len(panel);totals=counts.sum(axis=1,dtype=np.uint64);assert np.all(totals>0);logs=np.log1p(counts.astype('f8')/totals[:,None]*1e6)[:,idx];cohort=json.loads((src/'cohort.json').read_text());donors=cohort['eligibleDonors'];assert len(donors)==62
cohorts['GSE181897']={'donors':['GSE181897:exp_id:'+v for v in donors],'controls':logs[[lookup[(v,'C')] for v in donors]],'treated':logs[[lookup[(v,'B')] for v in donors]]}
arrays={}
for study,c in cohorts.items():
 for key in ['controls','treated']:arrays[study+'_'+key]=c[key]
np.savez_compressed(s/'cohorts.npz',**arrays);(s/'cohort-metadata.json').write_text(json.dumps({'featureIDs':panel,'studies':{study:c['donors'] for study,c in cohorts.items()}},indent=2)+'\n')
inputs=s/'inputs';inputs.mkdir(exist_ok=False)
for query in sorted(cohorts):
 train=[v for v in sorted(cohorts) if v!=query];donors=[d for v in train for d in cohorts[v]['donors']];studies=[v for v in train for _ in cohorts[v]['donors']];a={'featureIDs':panel,'trainingDonorIDs':donors,'trainingStudies':studies,'controls':np.concatenate([cohorts[v]['controls'] for v in train]).tolist(),'treated':np.concatenate([cohorts[v]['treated'] for v in train]).tolist(),'queryStudy':query,'queryDonorIDs':cohorts[query]['donors'],'queryControls':cohorts[query]['controls'].tolist()};(inputs/(query+'.json')).write_text(json.dumps(a,separators=(',',':'))+'\n')
(s/'input-freeze.json').write_text(json.dumps({'sourceBindings':bindings,'originalPredictionInputFreezeSHA256':sha(src/'input-freeze.json'),'inputs':{p.name:sha(p) for p in inputs.iterdir()},'cohortsSHA256':sha(s/'cohorts.npz'),'metadataSHA256':sha(s/'cohort-metadata.json'),'status':'three-study-inputs-prepared','sourceNativeGSECountsReconstructed':True,'panelFeatures':len(panel),'donors':{k:len(v['donors']) for k,v in cohorts.items()}},indent=2)+'\n')
print('Prepared', {k:len(v['donors']) for k,v in cohorts.items()},len(panel))
