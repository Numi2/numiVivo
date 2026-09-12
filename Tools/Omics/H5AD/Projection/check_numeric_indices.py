from pathlib import Path
import argparse,hashlib,json,subprocess,warnings,shutil
from importlib.metadata import version
import anndata as ad,h5py,numpy as np,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--legacy',type=Path,required=True);p.add_argument('--modern',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args();out=args.out;out.mkdir(exist_ok=False);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();cases=[]
def values(n,name):
 if name=='signed':return np.arange(n,dtype='i8')-2**62
 if name=='unsigned':return np.arange(n,dtype='u8')+np.uint64(2**63)
 if name=='big-endian':return (np.arange(n,dtype='i8')-10).astype('>i8')
 if name=='float':return np.array([0.,-0.,1.25,np.inf,-np.inf,np.nan,3.5][:n],dtype='f8')
 if name in ['boolean','nullable-boolean']:return np.arange(n)%2==0
 return np.arange(n,dtype='i8')-10
for mode in ['legacy','modern']:
 for name in ['signed','unsigned','big-endian','float','boolean']+(['nullable-integer','nullable-boolean'] if mode=='modern' else []):
  source=out/(mode+'-'+name+'.h5ad');shutil.copy2(args.legacy if mode=='legacy' else args.modern,source)
  with h5py.File(source,'r+') as f:
   for axis in ['obs','var','raw.var' if mode=='legacy' else 'raw/var']:
    if mode=='legacy':
     a=f[axis][()];names=a.dtype.names;v=values(len(a),name);b=np.empty(len(a),dtype=[(names[0],v.dtype)]+[(n,a.dtype.fields[n][0]) for n in names[1:]]);b[names[0]]=v
     for n in names[1:]:b[n]=a[n]
     del f[axis];f[axis]=b
    else:
     frame=f[axis];index=frame.attrs['_index'];n=len(frame[index]);del frame[index];v=values(n,name)
     if name.startswith('nullable-'):
      d=frame.create_group(index);d.attrs['encoding-type']=name;d.attrs['encoding-version']='0.1.0'
      for k,data in [('values',v),('mask',np.arange(n)==1)]:
       x=d.create_dataset(k,data=data);x.attrs['encoding-type']='array';x.attrs['encoding-version']='0.2.0'
     else:
      d=frame.create_dataset(index,data=v);d.attrs['encoding-type']='array';d.attrs['encoding-version']='0.2.0'
  with warnings.catch_warnings():warnings.simplefilter('ignore');original=ad.read_h5ad(source)
  for selection,rows,cols in [('full',None,None),('repeated',[6,0,3,0],[4,1,1,0]),('empty',[],[])]:
   tag=mode+'-'+name+'-'+selection;plan=out/(tag+'.json');dest=out/tag;plan.write_text(json.dumps({'schemaVersion':1,'source':{'bytes':list(bytes.fromhex(sha(source)))},'provenance':'Independent numeric index migration check; no biological claim','observationIndices':rows,'featureIndices':cols}))
   result=subprocess.run([str(args.binary),'project',str(source),str(plan),str(dest)],capture_output=True,text=True);(out/(tag+'.log')).write_text(result.stdout+result.stderr);assert result.returncode==0,(tag,result.stdout,result.stderr)
   with warnings.catch_warnings(),ad.settings.override(remove_unused_categories=False):
    warnings.simplefilter('ignore');expected=original[rows if rows is not None else slice(None),cols if cols is not None else slice(None)].copy();actual=ad.read_h5ad(dest/'projected.h5ad')
   for a,b in [(expected.obs,actual.obs),(expected.var,actual.var),(expected.raw.var,actual.raw.var)]:pd.testing.assert_frame_equal(a,b)
   for a,b in [(expected.X,actual.X),(expected.raw.X,actual.raw.X)]:np.testing.assert_equal(a.toarray() if hasattr(a,'toarray') else a,b.toarray() if hasattr(b,'toarray') else b)
   with h5py.File(source) as f,h5py.File(dest/'projected.h5ad') as g:
    for axis,target,indices in [('obs','obs',rows),('var','var',cols),('raw.var' if mode=='legacy' else 'raw/var','raw/var',None)]:
     after=g[target][g[target].attrs['_index']]
     if mode=='legacy':before=f[axis][()][f[axis].dtype.names[0]];pairs=[(before,after[()])]
     else:
      before=f[axis][f[axis].attrs['_index']];pairs=[(before[k][()],after[k][()]) for k in ['values','mask']] if isinstance(before,h5py.Group) else [(before[()],after[()])]
     for a,b in pairs:
      a=a[indices] if indices is not None else a;assert a.dtype==b.dtype and a.tobytes()==b.tobytes(),(tag,axis)
   assert sha(dest/'original.h5ad')==sha(source)
   replay=subprocess.run([str(args.binary),'verify-project',str(dest)],capture_output=True,text=True);(out/(tag+'-replay.log')).write_text(replay.stdout+replay.stderr);assert replay.returncode==0
   cases.append({'case':tag,'status':'passed','sourceSHA256':sha(source),'outputSHA256':sha(dest/'projected.h5ad'),'indexBytesExact':True,'replay':True})
report={'status':'passed','binarySHA256':sha(args.binary),'checkerSHA256':sha(Path(__file__)),'anndataVersion':version('anndata'),'cases':cases};(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print('Passed',len(cases),'cases')
