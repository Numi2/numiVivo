from pathlib import Path
import argparse,json,hashlib,subprocess,shutil,warnings
import h5py,numpy as np,anndata as ad,pandas as pd
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args();out=args.out;out.mkdir(exist_ok=False);cases=[];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def run(tag,*cmd,accept=True):
 x=subprocess.run([str(args.binary),*map(str,cmd)],capture_output=True,text=True);(out/(tag+'.log')).write_text(x.stdout+x.stderr);assert (x.returncode==0)==accept,(tag,x.stdout,x.stderr);return x
for tag in ['legacy-signed','legacy-unsigned','legacy-big-endian','legacy-float','legacy-boolean','modern-signed','modern-unsigned','modern-big-endian','modern-float','modern-boolean','modern-nullable-integer','modern-nullable-boolean']:
 source=out/(tag+'.h5ad');shutil.copy2(args.fixtures/(tag+'-full')/'projected.h5ad',source)
 with h5py.File(source,'r+') as f:
  # Construct bounded count fixtures; these are not measured biological data.
  for path in ['X/data','raw/X/data']:
   a=f[path][()];a[a>1000000]=17;f[path][...]=a
  for frame in ['obs','var']:
   x=f[frame][f[frame].attrs['_index']]
   if isinstance(x,h5py.Group) and 'mask' in x:x['mask'][...]=False
  frame=f['raw/var'];n=len(frame[frame.attrs['_index']])
  if isinstance(frame[frame.attrs['_index']],h5py.Group):n=len(frame[frame.attrs['_index']]['values'])
  for name in ['mapped_feature','mapped_name']:
   d=frame.create_dataset(name,data=np.array([f'gene-{i}' for i in range(n)],dtype=object),dtype=h5py.string_dtype());d.attrs['encoding-type']='string-array';d.attrs['encoding-version']='0.2.0'
  columns=list(frame.attrs['column-order'])+['mapped_feature','mapped_name'];del frame.attrs['column-order'];frame.attrs.create('column-order',np.array(columns,dtype=object),dtype=h5py.string_dtype())
 with warnings.catch_warnings():warnings.simplefilter('ignore');before=ad.read_h5ad(source)
 n,m=before.shape;edits=[]
 for path,values in [('obs/mapped_barcode',[f'cell-{i}' for i in range(n)]),('obs/mapped_sample',['sample']*n),('var/mapped_feature',[f'gene-{i}' for i in range(m)]),('var/mapped_name',[f'gene-{i}' for i in range(m)])]:edits.append({'path':path,'mode':'add','value':{'string':{'shape':[len(values)],'values':values}}})
 plan=out/(tag+'-annotation.json');plan.write_text(json.dumps({'schemaVersion':1,'source':{'bytes':list(bytes.fromhex(sha(source)))},'provenance':'Explicit synthetic identity columns; original numeric indices remain unchanged.','edits':edits}));annotated=out/(tag+'-annotated.h5ad');run(tag+'-annotate','annotate',source,plan,annotated)
 with warnings.catch_warnings():warnings.simplefilter('ignore');after=ad.read_h5ad(annotated)
 pd.testing.assert_frame_equal(before.obs,after.obs[before.obs.columns]);pd.testing.assert_frame_equal(before.var,after.var[before.var.columns]);pd.testing.assert_frame_equal(before.raw.var,after.raw.var)
 with h5py.File(source) as f,h5py.File(annotated) as g:
  datasets=[];f.visititems(lambda name,x:datasets.append(name) if isinstance(x,h5py.Dataset) else None)
  for name in datasets:assert f[name].dtype==g[name].dtype;np.testing.assert_equal(f[name][()],g[name][()])
 for matrix in ['X','raw/X']:
  label=tag+'-'+matrix.replace('/','-');mapping={'schemaVersion':1,'id':'explicit-identities','evidence':'synthetic','sourceDescription':'Format conformance fixture; no biological validation','countUnit':'umiCount','matrixPath':matrix,'sampleColumn':'mapped_sample','barcodeColumn':'mapped_barcode','featureIDColumn':'mapped_feature','featureNameColumn':'mapped_name','samples':[{'id':'sample','biologicalReplicateID':'replicate','condition':'control','batchID':'unreported','organism':'unreported'}]};mp=out/(label+'-mapping.json');mp.write_text(json.dumps(mapping));dest=out/label;run(label+'-count','count-store',annotated,mp,dest);run(label+'-replay','verify-count-store',dest)
  expected=(after.X if matrix=='X' else after.raw.X).tocsr();expected.sum_duplicates();expected.eliminate_zeros();coo=expected.tocoo();want=sorted(zip(coo.row.tolist(),coo.col.tolist(),coo.data.tolist()));got=np.fromfile(dest/'counts.bin',dtype=[('row','<u4'),('feature','<u4'),('count','<u8')]);assert sorted((int(x['row']),int(x['feature']),int(x['count'])) for x in got)==want
  metadata=json.loads((dest/'metadata.json').read_text());assert [c['barcode'] for c in metadata['cells']]==[f'cell-{i}' for i in range(n)];assert [g['id'] for g in metadata['features']]==[f'gene-{i}' for i in range(m)]
  ap=out/(label+'-aggregate.json');ap.write_text(json.dumps({'schemaVersion':1,'mapping':mapping,'contrasts':[]}));agg=out/(label+'-aggregate');run(label+'-aggregate','aggregate',annotated,ap,agg);run(label+'-aggregate-replay','verify-aggregate',agg)
  bulk=json.loads((agg/'report.json').read_text())['pseudobulk'];assert len(bulk['groups'])==1;bm=bulk['matrix'];dense=np.zeros(m,dtype='u8')
  for j,v in zip(bm['featureIndices'],bm['counts']):dense[j]=v
  np.testing.assert_equal(dense,np.asarray(expected.sum(axis=0)).ravel());cases.append({'case':label,'status':'passed','records':len(got),'countsSHA256':sha(dest/'counts.bin'),'annotatedSHA256':sha(annotated),'countReplay':True,'aggregateReplay':True})
(out/'summary.json').write_text(json.dumps({'status':'passed','cases':cases,'binarySHA256':sha(args.binary),'checkerSHA256':sha(Path(__file__))},indent=2)+'\n');print('Passed',len(cases),'complete analytical routes')
