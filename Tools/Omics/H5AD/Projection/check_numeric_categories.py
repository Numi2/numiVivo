from pathlib import Path
import hashlib,json,subprocess,warnings,shutil
from importlib.metadata import version
import anndata as ad,h5py,numpy as np,pandas as pd
import argparse
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out;out.mkdir(exist_ok=False);binary=a.binary;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();cases=[]
labels={'signed':np.array([-2**63,0,2**63-1],dtype='i8'),'unsigned':np.array([2**63+1,2**63+2,2**64-1],dtype='u8'),'big-endian':np.array([-3,2,7],dtype='>i4'),'float32':np.array([1.25,-0.,np.inf],dtype='f4'),'float64':np.array([-np.inf,-0.,3.5],dtype='f8'),'duplicate':np.array([1,1,3]),'signed-zero':np.array([0.,-0.,1.]),'nan':np.array([1.,np.nan,3.])}
# Pandas cannot construct a categorical index from a big-endian buffer.
# Normalize category byte order only in disposable reference copies; native output
# is checked separately for exact original dtype and bytes before replay.
def read_reference(path):
 copy=out/(path.parent.name+'-'+path.stem+'-reference.h5ad');shutil.copy2(path,copy)
 with h5py.File(copy,'r+') as h:
  names=[]
  h.visititems(lambda n,o: names.append(n) if isinstance(o,h5py.Dataset) and (n.endswith('_categories') or n.endswith('/categories')) and o.dtype.byteorder=='>' else None)
  for n in names:
   values=h[n][()];attrs=dict(h[n].attrs);del h[n];d=h.create_dataset(n,data=values.astype(values.dtype.newbyteorder('=')))
   for k,v in attrs.items():d.attrs[k]=v
 return ad.read_h5ad(copy)
for name,values in labels.items():
 source=out/(name+'.h5ad');shutil.copy2(a.source,source)
 with h5py.File(source,'r+') as f:del f['uns/donor_categories'];f['uns/donor_categories']=values
 rejected=None
 try:
  with warnings.catch_warnings():warnings.simplefilter('ignore');original=read_reference(source)
 except ValueError as e:rejected=str(e)
 for selection,rows,cols in [('full',None,None),('repeated',[6,0,3,0],[4,1,1,0]),('empty',[],[])]:
  tag=name+'-'+selection;plan=out/(tag+'.json');dest=out/tag;plan.write_text(json.dumps({'schemaVersion':1,'source':{'bytes':list(bytes.fromhex(sha(source)))},'provenance':'Independent numeric categorical migration check; no biological claim','observationIndices':rows,'featureIndices':cols}))
  result=subprocess.run([str(binary),'project',str(source),str(plan),str(dest)],capture_output=True,text=True);(out/(tag+'.log')).write_text(result.stdout+result.stderr)
  if rejected:
   expected='null legacy category labels' if name=='nan' else 'duplicate legacy category labels';assert result.returncode!=0 and expected in result.stdout+result.stderr and not dest.exists();cases.append({'case':tag,'status':'expected-rejection','referenceError':rejected});continue
  assert result.returncode==0,(tag,result.stdout,result.stderr)
  with ad.settings.override(remove_unused_categories=False):expected=original[rows if rows is not None else slice(None),cols if cols is not None else slice(None)].copy()
  actual=read_reference(dest/'projected.h5ad');pd.testing.assert_frame_equal(expected.obs,actual.obs);pd.testing.assert_frame_equal(expected.var,actual.var);pd.testing.assert_frame_equal(expected.raw.var,actual.raw.var)
  assert set(expected.uns)==set(actual.uns)
  for k in expected.uns:np.testing.assert_equal(expected.uns[k],actual.uns[k])
  with h5py.File(dest/'projected.h5ad') as h:
   for axis in ['obs','var']:
    got=h[axis+'/donor/categories'][()];assert got.dtype==values.dtype and got.tobytes()==values.tobytes(),(tag,axis)
  assert sha(dest/'original.h5ad')==sha(source)
  replay=subprocess.run([str(binary),'verify-project',str(dest)],capture_output=True,text=True);assert replay.returncode==0,(tag,replay.stdout,replay.stderr)
  cases.append({'case':tag,'status':'passed','sourceSHA256':sha(source),'outputSHA256':sha(dest/'projected.h5ad'),'categoryBytesExact':True,'replay':True,'referenceCategoryByteOrderNormalized':name=='big-endian'})
report={'status':'passed','binarySHA256':sha(binary),'checkerSHA256':sha(Path(__file__)),'anndataVersion':version('anndata'),'cases':cases};(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print('Passed',len(cases),'cases')
