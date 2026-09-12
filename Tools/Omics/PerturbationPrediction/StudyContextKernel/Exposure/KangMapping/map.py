from pathlib import Path
import json,hashlib
import anndata
r=Path('/Users/n/numivivo-kang-exposure-mapping-20260912');r.mkdir(exist_ok=False)
p=Path('/Users/n/numivivo-native-mnn-20260910/kang/pca');model=Path('/Users/n/numivivo-parse-context-evaluation-20260912/source-metadata.json')
a=json.loads((p/'metadata.json').read_text());m=json.loads(model.read_text());obs=anndata.read_h5ad(p/'original.h5ad',backed='r').obs
assert obs.index.tolist()==[x['barcode'] for x in a['cells']]
assert obs['native_sample'].astype(str).tolist()==[x['sampleID'] for x in a['cells']]
rows=[]
for donor in m['studies']['Kang']:
 source=donor.removeprefix('Kang:');samples=[x for x in a['samples'] if x['donorID']==source];assert len(samples)==2 and {x['condition'] for x in samples}=={'ctrl','stim'}
 record={'modelDonorID':donor,'sourceDonorID':source,'conditions':[]}
 for s in samples:
  cells=[x for x in a['cells'] if x['sampleID']==s['id']];b=sum(x.get('group')=='B cells' for x in cells);assert b>0
  record['conditions'].append({'condition':s['condition'],'nativeSampleID':s['id'],'fullCohortCells':len(cells),'sourceAnnotatedBCells':b,'batchID':s['batchID'],'primaryLibraryAccession':None})
 rows.append(record)
assert len(rows)==8 and sum(c['sourceAnnotatedBCells'] for x in rows for c in x['conditions'])==2651
h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
out={'status':'derived-source-identities-verified-primary-library-mapping-open','sources':{str(f):h(f) for f in [p/'metadata.json',p/'original.h5ad',model]},'rows':rows,'cellIdentityOrderVerified':True,'doseAssignedToModel':False,'limitations':['Native sample IDs are derived from the deposited benchmark annotations; they are not primary library accessions.','Batch remains unreported.','This does not reverify donor mean-expression vectors or establish cross-study participant independence.']}
(r/'mapping.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: eight model donors, sixteen derived condition IDs, 2651 source-annotated B cells; primary library mapping remains open')
