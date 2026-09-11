#!/usr/bin/env python3
"""Native legacy H5AD versus installed AnnData, including the complete original Kang cohort.
Format fixtures test interchange only, not biological prediction.
"""
import argparse,hashlib,json,os,subprocess,time,warnings
from pathlib import Path
from importlib.metadata import version
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--kang',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
summary=dict(status='in-progress',binarySHA256=sha(a.binary),checkerSHA256=sha(Path(__file__)),versions={k:version(k) for k in ['anndata','h5py','numpy','pandas','scipy']},cases=[],rejections=[])
def save():(a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
def equal(x,y):
 if sparse.issparse(x):
  assert sparse.issparse(y) and x.shape==y.shape and x.dtype==y.dtype
  x=x.tocsr(copy=True);y=y.tocsr(copy=True);x.sum_duplicates();y.sum_duplicates();x.sort_indices();y.sort_indices()
  for k in ('indptr','indices','data'):np.testing.assert_array_equal(getattr(x,k),getattr(y,k))
 elif isinstance(x,pd.DataFrame):pd.testing.assert_frame_equal(x,y)
 elif isinstance(x,dict):
  assert set(x)==set(y),(set(x),set(y))
  for k in x:equal(x[k],y[k])
 elif isinstance(x,np.ndarray):
  assert x.shape==y.shape and x.dtype==y.dtype,(x.shape,y.shape,x.dtype,y.dtype)
  np.testing.assert_array_equal(x,y)
  if x.dtype.kind not in 'OUS':assert np.array_equal(np.ascontiguousarray(x).view('u1'),np.ascontiguousarray(y).view('u1'))
 else:np.testing.assert_equal(x,y)
def run(name,source,rows=None,cols=None,reject=None):
 plan=dict(schemaVersion=1,source=dict(bytes=list(bytes.fromhex(sha(source)))),provenance='Complete legacy-source projection and AnnData comparison; no biological claim',observationIndices=rows,featureIndices=cols)
 pp=a.out/(name+'-plan.json');pp.write_text(json.dumps(plan)+'\n');dest=a.out/name
 with (a.out/(name+'.log')).open('w') as log:
  started=time.monotonic();r=subprocess.run([str(a.binary),'project',str(source),str(pp),str(dest)],stdout=log,stderr=subprocess.STDOUT)
  elapsed=time.monotonic()-started
  if reject:
   assert r.returncode>0 and not dest.exists(),(name,r.returncode)
   assert reject in (a.out/(name+'.log')).read_text(),name
   summary['rejections'].append(dict(name=name,sourceSHA256=sha(source),exitCode=r.returncode,diagnostic=reject));save();return
  assert r.returncode==0,(name,(a.out/(name+'.log')).read_text())
  subprocess.run([str(a.binary),'verify-project',str(dest)],stdout=log,stderr=subprocess.STDOUT,check=True)
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always');original=ad.read_h5ad(source)
  with ad.settings.override(remove_unused_categories=False):expected=original[rows if rows is not None else slice(None),cols if cols is not None else slice(None)].copy()
  actual=ad.read_h5ad(dest/'projected.h5ad')
 assert actual.shape==expected.shape
 for slot in ['X','obs','var','obsm','varm','obsp','varp','layers','uns']:
  x=getattr(expected,slot);y=getattr(actual,slot)
  equal(dict(x),dict(y)) if slot in ['obsm','varm','obsp','varp','layers','uns'] else equal(x,y)
 assert (expected.raw is None)==(actual.raw is None)
 if expected.raw is not None:
  equal(expected.raw.X,actual.raw.X);equal(expected.raw.var,actual.raw.var);equal(dict(expected.raw.varm),dict(actual.raw.varm))
 assert sha(dest/'original.h5ad')==sha(source)
 # Check original compound member representations as stored, including values
 # above float64 integer precision and fixed byte strings. Index/category
 # semantic decoding is separately compared through AnnData above.
 stored=0;sparse_sequences=0
 with h5py.File(source) as old,h5py.File(dest/'projected.h5ad') as new:
  for slot,selection in [('obs',rows),('var',cols),('raw.var',None),('obsm',rows),('varm',cols),('raw.varm',None)]:
   if slot not in old or not isinstance(old[slot],h5py.Dataset) or not old[slot].dtype.names:continue
   for field in old[slot].dtype.names:
    x=old[slot][field];x=x if selection is None else x[np.asarray(selection,dtype=int)]
    obj=new[slot.replace('raw.','raw/')+'/'+field];obj=obj['codes' if 'codes' in obj else 'values'] if isinstance(obj,h5py.Group) else obj
    equal(x,obj[()]);stored+=1
  for slot in ['X','raw.X']:
   if slot not in old or not isinstance(old[slot],h5py.Group):continue
   src=old[slot];dst=new[slot.replace('raw.','raw/')];kind=src.attrs['h5sparse_format'];dims=src.attrs['h5sparse_shape']
   r=np.arange(dims[0]) if rows is None else np.asarray(rows,dtype=int)
   c=np.arange(dims[1]) if cols is None or slot=='raw.X' else np.asarray(cols,dtype=int)
   major,minor=(r,c) if kind=='csr' else (c,r);minor_length=int(dims[1] if kind=='csr' else dims[0])
   multiplicity=np.bincount(minor,minlength=minor_length);offsets=np.r_[0,np.cumsum(multiplicity)];destinations=np.argsort(minor,kind='stable')
   ptr=src['indptr'][:];source_indices=src['indices'][:];source_data=src['data'][:];expected_ptr=[0];expected_indices=[];expected_data=[]
   for index in major:
    start,end=int(ptr[index]),int(ptr[index+1]);ii=source_indices[start:end];vv=source_data[start:end];copies=multiplicity[ii]
    positions=np.repeat(np.arange(len(ii)),copies);rank=np.arange(len(positions))-np.repeat(np.cumsum(copies)-copies,copies)
    mapped=destinations[offsets[ii[positions]]+rank];expected_indices.append(mapped);expected_data.append(vv[positions]);expected_ptr.append(expected_ptr[-1]+len(mapped))
   ii=np.concatenate(expected_indices) if expected_indices else np.array([],dtype='int64')
   vv=np.concatenate(expected_data) if expected_data else np.array([],dtype=src['data'].dtype)
   np.testing.assert_array_equal(expected_ptr,dst['indptr'][:]);np.testing.assert_array_equal(ii,dst['indices'][:]);equal(vv,dst['data'][:])
   assert dst['indices'].dtype==dst['indptr'].dtype==np.dtype('int64');sparse_sequences+=1
 summary['cases'].append(dict(name=name,sourceSHA256=sha(source),projectedSHA256=sha(dest/'projected.h5ad'),sourceShape=list(original.shape),outputShape=list(actual.shape),allAnnDataSlotsEqual=True,storedCompoundFieldsExact=stored,sparseStoredSequencesExact=sparse_sequences,reconstructionPassed=True,seconds=elapsed,warnings=sorted(set(str(w.message) for w in caught))))
 save();print(name,'PASS',actual.shape,flush=True)
def fixture(path,kind='csc',special=None):
 s=h5py.string_dtype('utf-8');n,m=7,5
 with h5py.File(path,'w') as f:
  obs=np.zeros(n,dtype=[('index',s),('donor','i1'),('huge','<u8'),('value','<f4'),('text','S8'),('skip','i1')])
  obs['index']=[f'c{i}' for i in range(n)];obs['donor']=[0,1,-1,0,1,0,1];obs['huge']=2**63+np.arange(n,dtype='u8');obs['value']=[1,np.nan,-np.inf,3,4,np.inf,-0.0];obs['text']=[b'a',b'b',b'c',b'd',b'e',b'f',b'g'];obs['skip']=[0,0,0,0,0,0,5]
  if special=='bad-code':obs['donor'][0]=-2
  if special in ['scalar','oversized-scalar']:obs['donor']=0
  f['obs']=obs;v=np.zeros(m,dtype=[('index',s),('name','S4'),('donor','i1')]);v['index']=[f'g{i}' for i in range(m)];v['name']=b'name';v['donor']=[0,1,0,1,0];f['var']=v
  x=np.arange(n*m,dtype='u8').reshape(n,m);x[x%3==0]=0;x[1,2]=2**63+17
  def matrix(name,data):
   if kind=='dense':f[name]=data;return
   z=(sparse.csr_matrix if kind=='csr' else sparse.csc_matrix)(data);g=f.create_group(name);g.attrs['h5sparse_format']=kind;g.attrs['h5sparse_shape']=z.shape
   for k in ['data','indices','indptr']:g[k]=getattr(z,k)
  matrix('X',x)
  e=np.zeros(n,dtype=[('pca','>f8',(3,)),('tensor','i4',(2,2))]);e['pca']=np.arange(n*3).reshape(n,3)/7;e['tensor']=np.arange(n*4).reshape(n,2,2);f['obsm']=e
  q=np.zeros(m,dtype=[('loadings','f4',(2,))]);q['loadings']=np.arange(m*2).reshape(m,2);f['varm']=q
  u=f.create_group('uns');u.create_dataset('donor_categories',data='a' if special=='scalar' else ['a','a'] if special=='duplicate' else ['b','a','unused'],dtype=s);u.create_dataset('skip_categories',data=['x'],dtype=s);u['retained']=np.arange(4)
  matrix('raw.X',x);f['raw.var']=v;f['raw.varm']=q
  if special=='mixed-raw':f.create_group('raw')
  if special=='numeric-categories':
   del u['donor_categories'];u['donor_categories']=np.array([1,2,3])
  if special=='oversized-scalar':
   del u['donor_categories'];u.create_dataset('donor_categories',data='x'*20000,dtype=s)
  if special=='frame-metadata':f['obs'].attrs['encoding-version']='0.2.0'
  if special=='scalar-embedding':
   del f['obsm'];f['obsm']=np.zeros(n,dtype=[('scalar','f8')])
  if special=='oversized-categories':
   del u['donor_categories'];u.create_dataset('donor_categories',shape=(100001,),dtype=s)

  if special=='partial':f.attrs['encoding-version']='0.1.0'
for kind in ['csr','csc','dense']:
 src=a.out/(kind+'.h5ad');fixture(src,kind)
 run(kind+'-full',src);run(kind+'-repeated',src,[6,0,3,0],[4,1,1,0]);run(kind+'-empty',src,[],[])
for special in ['scalar','bad-code','duplicate','mixed-raw','partial','numeric-categories','oversized-scalar','frame-metadata','scalar-embedding','oversized-categories']:
 src=a.out/(special+'.h5ad');fixture(src,special=special)
 reject={'bad-code':'below missing sentinel','duplicate':'duplicate legacy category labels','mixed-raw':'mixed legacy raw','partial':'encoding-type','numeric-categories':'requires string labels','oversized-scalar':'bounded legacy scalar string','frame-metadata':'contradictory legacy dataframe','scalar-embedding':'must be a fixed array','oversized-categories':'bounded scalar or vector'}.get(special)
 run(special,src,reject=reject)
assert sha(a.kang)=='e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830'
run('kang-original-complete',a.kang)
run('kang-original-repeated',a.kang,list(range(24672,-1,-1))+[0,24672],list(range(15705,-1,-1))+[0,15705])
summary['status']='passed';save()
