from pathlib import Path
import hashlib,json,subprocess,warnings,shutil,argparse
from importlib.metadata import version
import anndata as ad,h5py,numpy as np,pandas as pd
parser=argparse.ArgumentParser();parser.add_argument('--binary',type=Path,required=True);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();out=args.out;out.mkdir(exist_ok=False);binary=args.binary;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();cases=[]
for name,dtype in [('float','<f8'),('big-endian','>f4'),('signed','<i8'),('unsigned','<u8'),('boolean','?'),('string','S8')]:
 source=out/(name+'.h5ad');shutil.copy2(args.source,source)
 with h5py.File(source,'r+') as f:
  for key in ['obsm','varm','raw.varm']:
   n=len(f[key]);del f[key];a=np.zeros(n,dtype=[('scalar',dtype),('vector','<f4',(2,))])
   values=np.arange(n)
   if name=='unsigned':values=values.astype('u8')+np.uint64(2**63)
   elif name=='signed':values=values.astype('i8')-2**62
   elif name=='string':values=np.array([('cell-'+str(i)).encode() for i in range(n)],dtype=dtype)
   a['scalar']=values;a['vector']=np.arange(n*2).reshape(n,2);f[key]=a
 with warnings.catch_warnings():warnings.simplefilter('ignore');original=ad.read_h5ad(source)
 for selection,rows,cols in [('full',None,None),('repeated',[6,0,3,0],[4,1,1,0]),('empty',[],[])]:
  tag=name+'-'+selection;plan=out/(tag+'.json');dest=out/tag;plan.write_text(json.dumps({'schemaVersion':1,'source':{'bytes':list(bytes.fromhex(sha(source)))},'provenance':'Independent scalar embedding migration check; no biological claim','observationIndices':rows,'featureIndices':cols}))
  result=subprocess.run([str(binary),'project',str(source),str(plan),str(dest)],capture_output=True,text=True);(out/(tag+'.log')).write_text(result.stdout+result.stderr);assert result.returncode==0,(tag,result.stdout,result.stderr)
  with ad.settings.override(remove_unused_categories=False):expected=original[rows if rows is not None else slice(None),cols if cols is not None else slice(None)].copy()
  with warnings.catch_warnings():warnings.simplefilter('ignore');actual=ad.read_h5ad(dest/'projected.h5ad')
  pd.testing.assert_frame_equal(expected.obs,actual.obs);pd.testing.assert_frame_equal(expected.var,actual.var);pd.testing.assert_frame_equal(expected.raw.var,actual.raw.var)
  for a,b in [(expected.obsm,actual.obsm),(expected.varm,actual.varm),(expected.raw.varm,actual.raw.varm)]:
   assert set(a)==set(b)
   for k in a:np.testing.assert_equal(a[k],b[k])
  with h5py.File(source) as f,h5py.File(dest/'projected.h5ad') as g:
   for before,after,indices in [('obsm','obsm',rows),('varm','varm',cols),('raw.varm','raw/varm',None)]:
    values=f[before][()]['scalar'];values=values[indices] if indices is not None else values;got=g[after+'/scalar'][()];assert got.dtype==values.dtype and got.tobytes()==values.tobytes(),(tag,after)
  assert sha(dest/'original.h5ad')==sha(source)
  replay=subprocess.run([str(binary),'verify-project',str(dest)],capture_output=True,text=True);(out/(tag+'-replay.log')).write_text(replay.stdout+replay.stderr);assert replay.returncode==0,(tag,replay.stdout,replay.stderr)
  cases.append({'case':tag,'status':'passed','sourceSHA256':sha(source),'outputSHA256':sha(dest/'projected.h5ad'),'scalarBytesExact':True,'replay':True})
report={'status':'passed','binarySHA256':sha(binary),'checkerSHA256':sha(Path(__file__)),'anndataVersion':version('anndata'),'cases':cases};(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print('Passed',len(cases),'cases')
