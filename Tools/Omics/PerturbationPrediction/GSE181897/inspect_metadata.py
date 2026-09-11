from pathlib import Path
import sys,json,time,hashlib
r=Path(__file__).parent;sys.path.insert(0,str(r/'deps'))
import indexed_gzip,h5py,numpy as np
p=r/'GSE181897_concat.4.raw.h5ad.gz'; start=time.time()
assert p.is_file() and (r/'download.json').is_file()
with indexed_gzip.IndexedGzipFile(str(p),spacing=1024*1024) as gz:

 if (r/'source.gzidx').exists():gz.import_index(str(r/'source.gzidx'))
 else:gz.build_full_index();gz.export_index(str(r/'source.gzidx'))
 with h5py.File(gz,'r') as f:
  def conv(v):
   if isinstance(v,h5py.Reference):return {'referencedPath': f[v].name} if v else None
   if isinstance(v,bytes):return v.decode()
   if isinstance(v,np.ndarray):return [conv(x) for x in v.tolist()]
   if isinstance(v,np.generic):return conv(v.item())
   return v
  inventory=[]
  def visit(n,o):
   inventory.append(dict(path=n,type=type(o).__name__,attrs={k:conv(v) for k,v in o.attrs.items()},**(dict(shape=o.shape,dtype=str(o.dtype),chunks=o.chunks,compression=o.compression) if isinstance(o,h5py.Dataset) else {})))
  f.visititems(visit)
  (r/'h5-inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
  print(json.dumps(inventory,indent=2),flush=True)
print('elapsed',time.time()-start)
