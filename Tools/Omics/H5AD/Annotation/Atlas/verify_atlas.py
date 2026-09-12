from pathlib import Path
import argparse,json,hashlib,math,time,shutil
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--before',type=Path,required=True);p.add_argument('--report',type=Path,required=True);a=p.parse_args();before=json.loads(a.before.read_text());started=time.time();datasets=[];elements=0;groups=0
with h5py.File(a.source,'r') as f,h5py.File(a.output,'r') as g:
 names=[];f.visit(names.append)
 for name in ['']+names:
  x=f[name] if name else f;y=g[name] if name else g
  for key in x.attrs:
   if name=='obs' and key=='column-order':np.testing.assert_equal(list(y.attrs[key]),list(x.attrs[key])+['numivivo_source_row'])
   else:np.testing.assert_equal(x.attrs[key],y.attrs[key])
  if isinstance(x,h5py.Group):groups+=1;continue
  assert isinstance(y,h5py.Dataset) and x.shape==y.shape and x.dtype==y.dtype,name
  if x.shape==():pairs=[(x[()],y[()])]
  else:
   step=max(1,(8192 if x.dtype.hasobject else 262144)//max(1,math.prod(x.shape[1:])))
   pairs=((x[i:i+step],y[i:i+step]) for i in range(0,x.shape[0],step))
  for v,w in pairs:
   if x.dtype.hasobject:np.testing.assert_equal(v,w)
   else:assert np.asarray(v).tobytes()==np.asarray(w).tobytes(),name
  n=math.prod(x.shape);elements+=n;datasets.append({'path':name,'shape':x.shape,'dtype':str(x.dtype),'elements':n});print('checked',name,n,flush=True)
 row=g['obs/numivivo_source_row'];assert row.shape==(before['cells'],) and row.dtype==np.dtype('i8')
 assert row.attrs['encoding-type']=='array' and row.attrs['encoding-version']=='0.2.0'
 for i in range(0,before['cells'],65536):np.testing.assert_equal(row[i:i+65536],np.arange(i,min(i+65536,before['cells']),dtype='i8'))
 assert g['obs'].id!=g['var'].id
hashes={}
for name,path in [('source',a.source),('output',a.output)]:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 hashes[name]=h.hexdigest()
assert hashes['source']==before['sourceSHA256'];assert a.source.stat().st_ino!=a.output.stat().st_ino
report={'status':'passed','source':str(a.source),'output':str(a.output),'sourceSHA256':hashes['source'],'outputSHA256':hashes['output'],'sourceBytes':a.source.stat().st_size,'outputBytes':a.output.stat().st_size,'cellsAnnotated':before['cells'],'originalDatasetsChecked':datasets,'originalElementsChecked':elements,'originalGroupsChecked':groups,'sourceUnchanged':True,'separateInodes':True,'newColumnExact':True,'seconds':time.time()-started,'freeBytesAfter':shutil.disk_usage(a.output.parent).free,'scope':'Complete real HIRISA source preservation and row-provenance annotation; not biological labels or prediction validation.'};a.report.write_text(json.dumps(report,indent=2)+'\n');print('PASS',len(datasets),'datasets',elements,'elements',flush=True)
