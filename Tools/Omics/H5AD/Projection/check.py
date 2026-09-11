#!/usr/bin/env python3
"""Compare native axis projection with AnnData on format cases and real H5AD.

Format fixtures qualify interoperability only; no synthetic biological claim.
"""
import argparse, hashlib, json, os, shutil, subprocess, time
from pathlib import Path
from importlib.metadata import version
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--scoped',action='store_true')
p.add_argument('--out',type=Path,required=True);p.add_argument('--real',type=Path,action='append',default=[])
p.add_argument('--real-mapping',type=Path,action='append',default=[])
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
assert not a.real_mapping or len(a.real_mapping)==len(a.real)
summary=dict(status='in-progress',binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),cases=[],negativeCases=[],
 versions={k:version(k) for k in ['anndata','numpy','pandas','scipy','h5py']})
def save(): (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def invoke(source,plan,dest,label,verify=False,reject=False,existing=False):
 if a.scoped:args=['verify-project',str(dest)] if verify else ['project',str(source),str(plan),str(dest)]
 else:args=['singlecell-h5ad-project-verify',str(dest)] if verify else ['singlecell-h5ad-project',str(source),'--plan',str(plan),'--output',str(dest)]
 start=time.monotonic()
 with (a.out/(label+'.log')).open('w') as log:
  r=subprocess.run([str(a.binary.resolve())]+args,stdout=log,stderr=subprocess.STDOUT,
   env={**os.environ,'NUMIVIVO_HDF5_LIBRARY':os.environ.get('NUMIVIVO_HDF5_LIBRARY','/opt/homebrew/opt/hdf5/lib/libhdf5.dylib')})
 if reject:
  assert r.returncode!=0 and (verify or existing or not dest.exists()),(label,r.returncode)
 else: assert r.returncode==0,(label,r.returncode,(a.out/(label+'.log')).read_text()[-3000:])
 return dict(command=args,exitCode=r.returncode,seconds=time.monotonic()-start)
def plan_for(source,rows,cols):
 return dict(schemaVersion=1,source=dict(bytes=list(bytes.fromhex(digest(source)))),provenance='AnnData axis-projection interoperability comparison; original matrix semantics retained',
  observationIndices=rows,featureIndices=cols)
def equal(left,right):
 if sparse.issparse(left):
  assert sparse.issparse(right) and left.dtype==right.dtype and left.shape==right.shape
  x=left.tocsr(copy=True);y=right.tocsr(copy=True);x.sum_duplicates();y.sum_duplicates();x.sort_indices();y.sort_indices()
  for k in ['indptr','indices','data']:np.testing.assert_array_equal(getattr(x,k),getattr(y,k))
 elif isinstance(left,pd.DataFrame):pd.testing.assert_frame_equal(left,right)
 elif isinstance(left,dict):
  assert set(left)==set(right)
  for k in left:equal(left[k],right[k])
 elif isinstance(left,np.ndarray):
  assert isinstance(right,np.ndarray) and left.dtype==right.dtype and left.shape==right.shape
  np.testing.assert_array_equal(left,right)
  if left.dtype.kind not in 'OUS':assert np.array_equal(np.ascontiguousarray(left).view('u1'),np.ascontiguousarray(right).view('u1'))
 else:np.testing.assert_equal(left,right)
def compare(source,projected,rows,cols):
 original=ad.read_h5ad(source)
 with ad.settings.override(remove_unused_categories=False):expected=original[rows if rows is not None else slice(None),cols if cols is not None else slice(None)].copy()
 actual=ad.read_h5ad(projected)
 assert actual.shape==expected.shape
 for slot in ['X','obs','var','obsm','varm','obsp','varp','layers','uns']:
  left=getattr(expected,slot);right=getattr(actual,slot)
  equal(dict(left),dict(right)) if slot in ['obsm','varm','obsp','varp','layers','uns'] else equal(left,right)
 assert (expected.raw is None)==(actual.raw is None)
 if expected.raw is not None:
  assert actual.raw.shape==expected.raw.shape
  equal(expected.raw.X,actual.raw.X);equal(expected.raw.var,actual.raw.var);equal(dict(expected.raw.varm),dict(actual.raw.varm))
 # Storage-level dtype and attribute checks retain enum/string details and
 # nullable masked payloads that a high-level dataframe comparison may hide.
 storage=0
 with h5py.File(source) as src,h5py.File(projected) as dst:
  def axes(name):
   parts=name.split('/');categories='categories' in parts
   if parts[0] in ['obs','obsm'] and not categories:return rows,None
   if parts[0] in ['var','varm'] and not categories:return cols,None
   if parts[0]=='obsp':return rows,rows
   if parts[0]=='varp':return cols,cols
   if parts[0] in ['X','layers']:return rows,cols
   if parts[:2]==['raw','X']:return rows,None
   return None,None
  def visit(name,obj):
   nonlocal storage
   assert name in src,name
   old=src[name];assert type(old)==type(obj)
   for key in old.attrs:
    assert key in obj.attrs,(name,key)
    if key=='shape' and old.attrs.get('encoding-type') in ['csr_matrix','csc_matrix']:continue
    np.testing.assert_equal(old.attrs[key],obj.attrs[key])
   assert set(old.attrs)==set(obj.attrs)
   if not isinstance(obj,h5py.Dataset):
    kind=old.attrs.get('encoding-type')
    if kind in ['csr_matrix','csc_matrix'] and not name.startswith('uns/'):
     r,c=axes(name);dims=old.attrs['shape'];r=np.arange(dims[0]) if r is None else np.asarray(r,dtype=int);c=np.arange(dims[1]) if c is None else np.asarray(c,dtype=int)
     major,minor=(r,c) if kind=='csr_matrix' else (c,r);minor_length=dims[1] if kind=='csr_matrix' else dims[0]
     multiplicity=np.bincount(minor,minlength=minor_length);offsets=np.r_[0,np.cumsum(multiplicity)];destinations=np.argsort(minor,kind='stable')
     ptr=old['indptr'][:];expected_ptr=[0];expected_indices=[];expected_data=[]
     for row in major:
      start,end=int(ptr[row]),int(ptr[row+1]);index=old['indices'][start:end];values=old['data'][start:end];copies=multiplicity[index]
      source_positions=np.repeat(np.arange(len(index)),copies);rank=np.arange(len(source_positions))-np.repeat(np.cumsum(copies)-copies,copies)
      mapped=destinations[offsets[index[source_positions]]+rank]
      expected_indices.append(mapped);expected_data.append(values[source_positions]);expected_ptr.append(expected_ptr[-1]+len(mapped))
     ii=np.concatenate(expected_indices) if expected_indices else np.array([],dtype='int64')
     vv=np.concatenate(expected_data) if expected_data else np.array([],dtype=old['data'].dtype)
     np.testing.assert_array_equal(expected_ptr,obj['indptr'][:]);np.testing.assert_array_equal(ii,obj['indices'][:]);equal(vv,obj['data'][:])
    return
   is_sparse=old.parent.attrs.get('encoding-type') in ['csr_matrix','csc_matrix']
   if is_sparse and not name.startswith('uns/') and name.rsplit('/',1)[-1] in ['indices','indptr']:
    assert obj.dtype==np.dtype('int64');return
   assert old.dtype==obj.dtype,(name,old.dtype,obj.dtype)
   storage+=1
   if is_sparse:return  # Values compared as sparse matrices above.
   axis0,axis1=axes(name)
   values=old[()]
   if axis0 is not None:values=np.take(values,axis0,axis=0)
   if axis1 is not None:values=np.take(values,axis1,axis=1)
   equal(values,obj[()])
  dst.visititems(visit)
 return dict(sourceShape=list(original.shape),outputShape=list(actual.shape),rawShape=list(actual.raw.shape) if actual.raw is not None else None,
  storedDatasetDtypesChecked=storage,allAnnDataSlotsEqual=True)
def case(name,source,rows,cols,kind,**options):
 plan={**plan_for(source,rows,cols),**options};path=a.out/(name+'-plan.json');path.write_text(json.dumps(plan)+'\n');dest=a.out/name
 run=invoke(source,path,dest,name+'-publish')
 check=compare(source,dest/'projected.h5ad',rows,cols)
 verify=invoke(source,path,dest,name+'-verify',verify=True)
 summary['cases'].append(dict(name=name,kind=kind,sourceSHA256=digest(source),outputSHA256=digest(dest/'projected.h5ad'),run=run,verify=verify,**check));save()
 return source,plan

def fixture(path,kind):
 n,m=7,9
 x=np.arange(n*m,dtype='uint64').reshape(n,m);x[x%3==0]=0;x[1,2]=2**53+3;x[5,8]=2**63+7
 if kind=='csr':matrix=sparse.csr_matrix(x)
 elif kind=='csc':matrix=sparse.csc_matrix(x.astype('float64')/7)
 else:matrix=(x%37).astype('>f4');matrix[1,2]=np.nan;matrix[4,3]=np.inf
 obs=pd.DataFrame(index=pd.Index([f'cell-{i}' for i in range(n)],dtype='string'))
 obs['donor']=pd.Categorical(['d1','d1','d2',None,'d3','d3','d4'],categories=['d4','d3','d2','d1','unused'],ordered=True)
 obs['integer']=pd.array([2**60,None,-2**60,3,4,5,6],dtype='Int64')
 obs['boolean']=pd.array([True,None,False,True,False,True,None],dtype='boolean')
 obs['text']=pd.array(['α',None,'β','🙂','x','long',''],dtype='string')
 obs['float']=np.array([1,np.nan,-np.inf,0,3,4,np.inf],dtype='float32')
 var=pd.DataFrame({'symbol':pd.Categorical([f'g{i%3}' for i in range(m)])},index=[f'gene-{i}' for i in range(m)])
 data=ad.AnnData(matrix,obs=obs,var=var)
 data.layers['fractional']=sparse.csc_matrix(np.arange(n*m,dtype='float32').reshape(n,m)/11)
 data.layers['dense']=np.arange(n*m,dtype='int16').reshape(n,m)
 data.obsm['pca']=np.arange(n*4,dtype='float32').reshape(n,4)
 data.obsm['tensor']=np.arange(n*6,dtype='int32').reshape(n,2,3)
 data.obsm['frame']=pd.DataFrame({'nullable':pd.array([1,None,3,4,5,6,7],dtype='Int64'),'label':obs.donor.to_numpy()},index=obs.index)
 data.varm['loadings']=np.arange(m*3,dtype='float64').reshape(m,3)
 data.obsp['dense']=np.eye(n,dtype='bool');data.obsp['csc']=sparse.csc_matrix(np.eye(n,dtype='float32'))
 data.varp['complex']=np.eye(m,dtype='complex64')*(1+2j)
 raw=ad.AnnData(sparse.csr_matrix(np.arange(n*12,dtype='int32').reshape(n,12)),obs=pd.DataFrame(index=obs.index),var=pd.DataFrame(index=[f'raw-{i}' for i in range(12)]))
 raw.varm['raw-loadings']=np.arange(24,dtype='float32').reshape(12,2);data.raw=raw
 data.uns['nested']={'array':np.array([2**63+9],dtype='uint64'),'strings':np.array(['a','β']),'scalar':np.float32(np.nan)}
 with ad.settings.override(allow_write_nullable_strings=True):data.write_h5ad(path,compression='gzip')
 with h5py.File(path,'r+') as f:
  f['uns']['original_X_alias']=f['X']
 return path

for kind in ['csr','csc','dense']:
 source=fixture(a.out/(kind+'-source.h5ad'),kind)
 case(kind+'-reorder',source,[5,1,6,0],[8,2,0,6], 'format-regression')
 case(kind+'-repeat',source,[5,1,5,6,1,0],[8,2,8,0,2,6], 'repeated-axis-format')
case('empty-cells',source,[],[3,1],'format-regression')
case('empty-features',source,[4,0],[],'format-regression')
case('identity',source,None,None,'format-regression')
empty=ad.AnnData(np.zeros((2,3),dtype='float32'));empty.obsm['embedding']=np.ones((2,4),dtype='float32')
empty_source=a.out/'empty-budget-source.h5ad';empty.write_h5ad(empty_source)
case('empty-budget',empty_source,[],[],'format-regression',maximumElementVisits=1)
# Unsorted duplicate sparse entries remain data, not implicit canonicalization.
source=fixture(a.out/'duplicates-source.h5ad','csr')
with h5py.File(source,'r+') as f:
 indices=f['X/indices'];values=indices[:];values[:3]=[2,1,2];indices[:]=values
case('sparse-duplicates',source,[0,5,1],[2,8,1,0],'format-regression')
case('sparse-duplicates-repeat',source,[0,5,0,1],[2,8,2,1,0,2],'repeated-axis-format')
# A single sparse source entry expands beyond the 65,536-value transfer bound.
for kind in ('csr','csc','dense'):
 values=np.array([[11,0,7],[0,4,0]],dtype='uint64')
 values[0,0]=2**53+3
 matrix=sparse.csr_matrix(values) if kind=='csr' else sparse.csc_matrix(values) if kind=='csc' else values
 large=a.out/(kind+'-repeat-buffer-source.h5ad');ad.AnnData(matrix).write_h5ad(large)
 rr,cc=([1,0,1],[1]*65537+[2,0,1]) if kind!='csc' else ([1]*65537+[0,1,0],[1,0,1])
 case(kind+'-repeat-buffer',large,rr,cc,'repeated-axis-transfer-boundary')
for i,source in enumerate(a.real):
 with h5py.File(source) as f:
  def length(frame):
   index=frame[frame.attrs['_index']]
   if isinstance(index,h5py.Dataset):return len(index)
   return len(index['codes'] if 'codes' in index else index['values'])
  n=length(f['obs']);m=length(f['var'])
 rows=list(range(n-1,-1,-2));cols=list(range(m-1,-1,-3))
 case('real-'+str(i),source,rows,cols,'real-data-interoperability')
 if a.real_mapping:
  mapping=json.loads(a.real_mapping[i].read_text());mapping=mapping.get('mapping',mapping)
  mapping_path=a.out/('real-'+str(i)+'-mapping.json');mapping_path.write_text(json.dumps(mapping)+'\n')
  projected=a.out/('real-'+str(i))/'projected.h5ad';imported=a.out/('real-'+str(i)+'-import')
  args=[str(projected),str(mapping_path),str(imported)] if a.scoped else ['singlecell-h5ad-import',str(projected),'--plan',str(mapping_path),'--output',str(imported)]
  with (a.out/('real-'+str(i)+'-import.log')).open('w') as log:
   result=subprocess.run([str(a.binary.resolve())]+args,stdout=log,stderr=subprocess.STDOUT)
  assert result.returncode==0,('count reimport',i,result.returncode)
  dataset=json.loads((imported/'dataset.json').read_text());expected=ad.read_h5ad(projected)
  assert [c['barcode'] for c in dataset['cells']]==list(expected.obs_names)
  assert [c['sampleID'] for c in dataset['cells']]==list(expected.obs[mapping['sampleColumn']].astype(str))
  assert [f['id'] for f in dataset['features']]==list(expected.var_names)
  mtx=dataset['matrix'];native=sparse.csr_matrix((np.array(mtx['counts'],dtype='uint64'),mtx['featureIndices'],mtx['rowOffsets']),shape=expected.shape)
  original=expected.X.tocsr();assert np.all(original.data==np.floor(original.data))
  equal(original.astype('uint64'),native)
  assert digest(imported/'original.h5ad')==digest(projected)
  summary['cases'][-1]['nativeCountReimport']=dict(exactCounts=True,exactAxes=True,exactOriginalH5AD=True,nonzeros=native.nnz,command=args)
  save()
# Retain every original real cell/feature, reorder and repeat declared positions.
for i,source in enumerate(a.real):
 original=ad.read_h5ad(source);n,m=original.shape
 rr=list(range(n-1,-1,-1))+[0,n//2,n-1];cc=list(range(m-1,-1,-1))+[0,m//2,m-1]
 case('real-repeat-'+str(i),source,rr,cc,'complete-source-repeated-axis-interoperability')
 del original
# Rejection must never publish a partial destination.
base=plan_for(source,[0],[0])
for name,edit in [('negative-index',{'observationIndices':[-1,0]}),('out-of-bounds',{'featureIndices':[1999999]}),('work-limit',{'maximumElementVisits':1}),('zero-storage',{'maximumOutputBytes':0}),('negative-storage',{'maximumOutputBytes':-1}),('excess-storage',{'maximumOutputBytes':8589934593}),('tiny-storage',{'maximumOutputBytes':1}),('wrong-source',{'source':{'bytes':[0]*32}})]:
 plan={**base,**edit};path=a.out/(name+'-plan.json');path.write_text(json.dumps(plan)+'\n')
 summary['negativeCases'].append(dict(name=name,**invoke(source,path,a.out/name,name,reject=True)));save()
expanded=a.out/'expanded-source.h5ad';ad.AnnData(sparse.csr_matrix([[1]],dtype='uint64')).write_h5ad(expanded)
expanded_plan=a.out/'expanded-plan.json';expanded_plan.write_text(json.dumps(plan_for(expanded,[0]*70000,[0]*70000))+'\n')
summary['negativeCases'].append(dict(name='repeated-sparse-work-limit',**invoke(expanded,expanded_plan,a.out/'expanded-output','expanded-work',reject=True)))
# Unknown aligned encoding cannot be silently retained with the original axis.
bad=a.out/'unknown-encoding.h5ad';shutil.copy2(a.out/'csr-source.h5ad',bad)
with h5py.File(bad,'r+') as f:f['obsm/pca'].attrs['encoding-type']='unknown'
path=a.out/'unknown-plan.json';path.write_text(json.dumps(plan_for(bad,[0],[0]))+'\n')
summary['negativeCases'].append(dict(name='unknown-encoding',**invoke(bad,path,a.out/'unknown','unknown',reject=True)))
# A malformed aligned matrix must not publish a partially transformed file.
bad_shape=a.out/'misaligned.h5ad';shutil.copy2(a.out/'csr-source.h5ad',bad_shape)
with h5py.File(bad_shape,'r+') as f:
 del f['layers/dense'];d=f['layers'].create_dataset('dense',data=np.zeros((3,3),dtype='float32'))
 d.attrs['encoding-type']='array';d.attrs['encoding-version']='0.2.0'
bad_plan=a.out/'misaligned-plan.json';bad_plan.write_text(json.dumps(plan_for(bad_shape,[0],[0]))+'\n')
summary['negativeCases'].append(dict(name='misaligned',**invoke(bad_shape,bad_plan,a.out/'misaligned-output','misaligned',reject=True)))
bad_column=a.out/'matrix-column.h5ad';shutil.copy2(a.out/'csr-source.h5ad',bad_column)
with h5py.File(bad_column,'r+') as f:
 del f['obs/float'];d=f['obs'].create_dataset('float',data=np.zeros((7,2),dtype='float32'))
 d.attrs['encoding-type']='array';d.attrs['encoding-version']='0.2.0'
column_plan=a.out/'matrix-column-plan.json';column_plan.write_text(json.dumps(plan_for(bad_column,[0],[0]))+'\n')
summary['negativeCases'].append(dict(name='matrix-dataframe-column',**invoke(bad_column,column_plan,a.out/'matrix-column-output','matrix-column',reject=True)))
# A valid lazily allocated dense matrix can be tiny on disk but exceed the
# output allowance. Reject before allocating/writing its 1.6 GB value array.
oversized=a.out/'oversized-dense.h5ad';ad.AnnData(shape=(10000,20000)).write_h5ad(oversized)
with h5py.File(oversized,'r+') as f:
 d=f.create_dataset('X',shape=(10000,20000),dtype='float64',chunks=(100,100),compression='gzip')
 d.attrs['encoding-type']='array';d.attrs['encoding-version']='0.2.0'
oversized_plan=a.out/'oversized-plan.json';oversized_plan.write_text(json.dumps(plan_for(oversized,None,None))+'\n')
summary['negativeCases'].append(dict(name='output-allocation-limit',**invoke(oversized,oversized_plan,a.out/'oversized-output','oversized',reject=True)))
# Explicit byte allowance changes admission/receipt, not projection values.
explicit_source=a.out/'csr-source.h5ad';explicit_plan=json.loads((a.out/'csr-reorder-plan.json').read_text());explicit_plan['maximumOutputBytes']=4*1024*1024
explicit_path=a.out/'explicit-storage-plan.json';explicit_path.write_text(json.dumps(explicit_plan)+'\n');explicit_out=a.out/'explicit-storage'
record=invoke(explicit_source,explicit_path,explicit_out,'explicit-storage')
assert digest(explicit_out/'projected.h5ad')==digest(a.out/'csr-reorder/projected.h5ad')
assert json.loads((explicit_out/'plan.json').read_text())['maximumOutputBytes']==4*1024*1024
assert 'maximumOutputBytes' not in json.loads((a.out/'csr-reorder/plan.json').read_text())
record['verify']=invoke(explicit_source,explicit_path,explicit_out,'explicit-storage-verify',verify=True)
summary['cases'].append(dict(name='explicit-storage',defaultOutputBytesExact=True,**record));save()
# Explicit no-overwrite checks the existing valid output remains unchanged.
existing=a.out/'csr-reorder';before=digest(existing/'projected.h5ad')
summary['negativeCases'].append(dict(name='no-overwrite',**invoke(a.out/'csr-source.h5ad',a.out/'csr-reorder-plan.json',existing,'no-overwrite',reject=True,existing=True)))
assert digest(existing/'projected.h5ad')==before
# Output tampering and no-overwrite behavior are exercised on a copied bundle.
bundle=a.out/'tampered';shutil.copytree(a.out/'csr-reorder',bundle)
with h5py.File(bundle/'projected.h5ad','r+') as f:f['X/data'][0]=17
summary['negativeCases'].append(dict(name='tampered-output',**invoke(source,path,bundle,'tampered',verify=True,reject=True)))
summary['status']='passed-H5AD-projection-interoperability';summary['qualification']='Exact interoperability on declared real files plus storage fixtures; no biological or million-cell qualification.'
save();print(json.dumps(summary,indent=2))
