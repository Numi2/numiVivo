#!/usr/bin/env python3
"""Qualify complete original Kang projection -> annotation -> counts/pseudobulk.

Python prepares explicit plans and checks native artifacts. It never writes a
replacement H5AD or a reduced cohort. This is data/analysis input qualification,
not a new model fit or independent biological prediction result.
"""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
from importlib.metadata import version
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--mapping',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def read(path):return json.loads(path.read_text())
def identity(path):return dict(bytes=list(bytes.fromhex(sha(path))))
summary=dict(status='in-progress',sourceSHA256=sha(a.source),binarySHA256=sha(a.binary),checkerSHA256=sha(Path(__file__)),mappingSHA256=sha(a.mapping),versions={k:version(k) for k in ['anndata','h5py','numpy','pandas','scipy']},commands=[])
def save():write(a.out/'summary.json',summary)
def run(name,*args):
 start=time.monotonic()
 with (a.out/(name+'.stdout')).open('wb') as out,(a.out/(name+'.log')).open('wb') as log:r=subprocess.run([str(a.binary),*map(str,args)],stdout=out,stderr=log)
 record=dict(name=name,arguments=list(map(str,args)),returncode=r.returncode,seconds=time.monotonic()-start);summary['commands'].append(record);save()
 assert r.returncode==0,(name,r.returncode,(a.out/(name+'.log')).read_text()[-2000:])
 return a.out/(name+'.stdout')
assert summary['sourceSHA256']=='e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830'
original=ad.read_h5ad(a.source);assert original.shape==(24673,15706)
raw=original.X.data;assert np.all(np.isfinite(raw)&(raw>=0)&(raw==np.floor(raw)))
counts=original.X.astype(np.uint64).tocsc();counts.sum_duplicates();counts.eliminate_zeros();counts.sort_indices()
assert counts.nnz==14184532
row_totals=np.asarray(counts.sum(axis=1)).ravel();column_totals=np.asarray(counts.sum(axis=0)).ravel();row_nnz=np.bincount(counts.indices,minlength=counts.shape[0]);column_nnz=np.diff(counts.indptr)
np.testing.assert_array_equal(row_totals,original.obs['nCount_RNA']);np.testing.assert_array_equal(row_nnz,original.obs['nFeature_RNA'])
mapping=read(a.mapping);mapping=mapping.get('mapping',mapping)
assert mapping['sampleColumn']=='native_sample' and mapping['groupColumn']=='cell_type' and mapping['matrixPath']=='X' and mapping['countUnit']=='umiCount'
roles={s['id']:s for s in mapping['samples']};sample=original.obs['replicate'].astype(str)+'__'+original.obs['label'].astype(str)
assert set(sample)==set(roles) and len(roles)==16
for donor,condition,sid in zip(original.obs['replicate'].astype(str),original.obs['label'].astype(str),sample):
 role=roles[sid];assert role['donorID']==role['biologicalReplicateID']==donor and role['condition']==condition and role['batchID']=='unreported'
write(a.out/'mapping.json',mapping)
plan=dict(schemaVersion=1,source=identity(a.source),provenance='Complete original Kang source preserved for explicitly mapped RNA-count analysis; all original cells, genes and annotations retained')
write(a.out/'projection-plan.json',plan)
run('projection','project',a.source,a.out/'projection-plan.json',a.out/'projection');run('projection-replay','verify-project',a.out/'projection')
projected=a.out/'projection/projected.h5ad'
with h5py.File(projected) as f:assert dict(f['uns'].attrs)=={'encoding-type':'dict','encoding-version':'0.1.0'}
annotation=dict(schemaVersion=1,source=identity(projected),provenance='Join each original author replicate and label with __ as an explicit sample identity; batch remains unreported. Preserve all original fields and category definitions.',edits=[dict(path='obs/native_sample',mode='add',value={'string':dict(shape=[len(sample)],values=sample.tolist())})])
write(a.out/'annotation-plan.json',annotation)
annotated=a.out/'annotated.h5ad';receipt=run('annotation','annotate',projected,a.out/'annotation-plan.json',annotated)
r=read(receipt);assert r['source']==identity(projected) and r['output']==identity(annotated)
def equal(x,y):
 if sparse.issparse(x):
  assert sparse.issparse(y) and x.shape==y.shape and x.dtype==y.dtype
  x=x.tocsr(copy=True);y=y.tocsr(copy=True);x.sum_duplicates();y.sum_duplicates();x.sort_indices();y.sort_indices()
  for k in ['data','indices','indptr']:np.testing.assert_array_equal(getattr(x,k),getattr(y,k))
 elif isinstance(x,pd.DataFrame):pd.testing.assert_frame_equal(x,y)
 elif isinstance(x,dict):
  assert set(x)==set(y)
  for k in x:equal(x[k],y[k])
 else:np.testing.assert_equal(x,y)
actual=ad.read_h5ad(annotated);np.testing.assert_array_equal(actual.obs.pop('native_sample'),sample);journal=actual.uns.pop('numivivo_edits');assert len(journal)==1
for entry in journal.values():
 assert entry['source_sha256']==sha(projected) and json.loads(entry['plan_json'])==annotation
for slot in ['X','obs','var','obsm','varm','obsp','varp','layers','uns']:
 x=getattr(original,slot);y=getattr(actual,slot);equal(dict(x),dict(y)) if slot in ['obsm','varm','obsp','varp','layers','uns'] else equal(x,y)
assert actual.raw is None and original.raw is None
stored_datasets=0
with h5py.File(projected) as src,h5py.File(annotated) as dst:
 def check_storage(name,obj):
  global stored_datasets
  target=dst[name];assert type(target)==type(obj)
  for key,value in obj.attrs.items():
   if name=='obs' and key=='column-order':np.testing.assert_array_equal(target.attrs[key],list(value)+['native_sample'])
   else:np.testing.assert_equal(value,target.attrs[key])
  if isinstance(obj,h5py.Dataset):
   assert obj.dtype==target.dtype and obj.shape==target.shape;np.testing.assert_array_equal(obj[()],target[()]);stored_datasets+=1
 src.visititems(check_storage)
run('count-store','count-store',annotated,a.out/'mapping.json',a.out/'count-store');run('count-store-replay','verify-count-store',a.out/'count-store')
store=a.out/'count-store';metadata=read(store/'metadata.json');quality=read(store/'quality.json');receipt=read(store/'receipt.json')
assert metadata['samples']==mapping['samples'];assert [f['id'] for f in metadata['features']]==original.var_names.tolist()
for i,cell in enumerate(metadata['cells']):assert cell==dict(barcode=str(original.obs_names[i]),sampleID=sample.iloc[i],group=str(original.obs['cell_type'].iloc[i]))
assert len(metadata['cells'])==counts.shape[0] and len(metadata['features'])==counts.shape[1]
for key,values in [('rowTotals',row_totals),('rowNonzeros',row_nnz),('featureTotals',column_totals),('featureNonzeros',column_nnz)]:np.testing.assert_array_equal(quality[key],values)
assert receipt['entries']==counts.nnz and receipt['source']==identity(annotated) and sha(store/'original.h5ad')==sha(annotated)
records=np.memmap(store/'counts.bin',mode='r',dtype=[('row','<u4'),('feature','<u4'),('count','<u8')]);assert len(records)==counts.nnz
for feature in range(counts.shape[1]):
 start,end=map(int,counts.indptr[feature:feature+2]);chunk=records[start:end]
 np.testing.assert_array_equal(chunk['row'],counts.indices[start:end]);np.testing.assert_array_equal(chunk['count'],counts.data[start:end]);assert np.all(chunk['feature']==feature)
del records
write(a.out/'aggregate-plan.json',dict(schemaVersion=1,mapping=mapping,contrasts=[]))
run('aggregate','aggregate',annotated,a.out/'aggregate-plan.json',a.out/'aggregate');run('aggregate-replay','verify-aggregate',a.out/'aggregate')
report=read(a.out/'aggregate/report.json');bulk=report['pseudobulk'];assert report['metadata']==metadata and report['canonicalNonzeros']==counts.nnz
keys=list(zip(original.obs['replicate'].astype(str),original.obs['label'].astype(str),original.obs['cell_type'].astype(str)));groups=sorted(set(keys));lookup={k:i for i,k in enumerate(groups)}
assignment=np.array([lookup[k] for k in keys]);incidence=sparse.csr_matrix((np.ones(len(keys),dtype=np.uint64),(assignment,np.arange(len(keys)))),shape=(len(groups),len(keys)))
expected=(incidence@counts).tocsr();expected.sum_duplicates();expected.sort_indices()
for i,(group,key) in enumerate(zip(bulk['groups'],groups)):
 members=np.flatnonzero(assignment==i).tolist();assert group['sourceCellIndices']==members
 assert group['biologicalReplicateID']==group['donorID']==key[0] and group['condition']==key[1] and group['cellGroup']==key[2]
 assert group['sampleIDs']==[key[0]+'__'+key[1]] and group['batchIDs']==['unreported']
assert len(bulk['groups'])==len(groups) and bulk['featureIDs']==original.var_names.tolist()
for key,values in [('rowOffsets',expected.indptr),('featureIndices',expected.indices),('counts',expected.data)]:np.testing.assert_array_equal(bulk['matrix'][key],values)
summary.update(status='passed',cells=counts.shape[0],features=counts.shape[1],nativeCountRecords=counts.nnz,samples=len(roles),donors=original.obs['replicate'].nunique(),sourceCellTypes=original.obs['cell_type'].nunique(),pseudobulkGroups=len(groups),pseudobulkNonzeros=expected.nnz,allCountRecordsExact=True,allSourceAnnotationsPreserved=True,storedDatasetsPreserved=stored_datasets,allCellAndFeatureIdentitiesExact=True,allDonorTreatmentAssignmentsExact=True,allPseudobulkValuesExact=True,allRNAAnnotationTotalsAndFeaturesExact=True,projectionSHA256=sha(projected),annotatedSHA256=sha(annotated),countStoreSHA256=sha(store/'counts.bin'),scope='Complete author-annotated Kang input and exact count/pseudobulk semantics. UMI designation follows the frozen earlier source mapping; RNA totals/nonzeros independently match every source cell. No new held-out prediction, DE model fit or general biological validation.')
save();print(json.dumps(summary,indent=2))
