from pathlib import Path
import argparse,h5py,numpy as np,json,hashlib,subprocess
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False);rows=[]
for n in [2000000,2000001]:
 source=a.out/(str(n)+'.h5ad')
 with h5py.File(source,'w') as f:
  f.attrs['encoding-type']='anndata';f.attrs['encoding-version']='0.1.0'
  for key,size in [('obs',n),('var',1)]:
   g=f.create_group(key);g.attrs['encoding-type']='dataframe';g.attrs['encoding-version']='0.2.0';g.attrs['_index']='index';g.attrs.create('column-order',np.array([],dtype=object),dtype=h5py.string_dtype());d=g.create_dataset('index',shape=(size,),dtype='i8');d.attrs['encoding-type']='array';d.attrs['encoding-version']='0.2.0'
  g=f.create_group('X');g.attrs['encoding-type']='csr_matrix';g.attrs['encoding-version']='0.1.0';g.attrs['shape']=np.array([n,1],dtype='i8');g.create_dataset('indptr',shape=(n+1,),dtype='i8');g.create_dataset('indices',shape=(0,),dtype='i8');g.create_dataset('data',shape=(0,),dtype='i8')
 plan=a.out/(str(n)+'.json');plan.write_text(json.dumps({'schemaVersion':1,'source':{'bytes':list(hashlib.sha256(source.read_bytes()).digest())},'provenance':'Software-only annotation axis boundary; unallocated zero-filled count fixture','edits':[{'path':'uns/check','mode':'add','value':{'string':{'shape':[],'values':['axis boundary']}}}]}));dest=a.out/(str(n)+'-annotated.h5ad');r=subprocess.run([str(a.binary),'annotate',str(source),str(plan),str(dest)],capture_output=True,text=True);(a.out/(str(n)+'.log')).write_text(r.stdout+r.stderr)
 if n==2000000:assert r.returncode==0 and dest.exists()
 else:assert r.returncode==65 and 'annotation dataframe index shape' in r.stderr and not dest.exists()
 rows.append({'rows':n,'status':'accepted' if r.returncode==0 else 'expected-rejection','sourceBytes':source.stat().st_size})
(a.out/'summary.json').write_text(json.dumps({'status':'passed','cases':rows},indent=2)+'\n');print(rows)
