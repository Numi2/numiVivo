#!/usr/bin/env python3
"""Source-bound admission into the existing native count store/composite corpus.

`freeze` reads metadata only. `admit` preserves sparse integer counts and invokes
native owners. A full source may be discarded after its digest and selected rows
are verified; the receipt distinguishes the original source from the projection.
"""
import argparse, collections, hashlib, json, subprocess
from pathlib import Path
import h5py
import numpy as np
from scipy import sparse
from wetlab import read, write, sha, require

TASKS=('known-target','held-target','held-context')
def fingerprint(hexadecimal):return {'bytes':list(bytes.fromhex(hexadecimal))}
def digest(value):return hashlib.sha256(value.encode()).hexdigest()
def column(group,name):
 d=group[name]
 if isinstance(d,h5py.Group):
  values=d['categories'][:];codes=d['codes'][:]
  return [None if n<0 else text(values[n]) for n in codes]
 values=d[:]
 if 'categories' in d.attrs:
  cats=group.file[d.attrs['categories']][:]
  return [None if int(n)<0 else text(cats[n]) for n in values]
 return [text(x) for x in values]
def text(x):return x.decode() if isinstance(x,bytes) else x.item() if hasattr(x,'item') else x

def fetch(manifest,source_id,output):
 """Bounded public-file transfer; whole-source digest verified before admission."""
 import concurrent.futures,os,requests,shutil,time
 spec=next(x for x in read(manifest)['sources'] if x['id']==source_id)
 out=Path(output);require(not out.exists(),'Download destination exists')
 url=spec.get('downloadURL',spec['url'])
 response=requests.get(url,headers={'Range':'bytes=0-0'},timeout=60)
 response.raise_for_status();require(response.status_code==206,'Source must support bounded range reads')
 size=int(response.headers['Content-Range'].split('/')[-1]);require(shutil.disk_usage(out.parent).free>size+1024**3,'Insufficient staging disk; existing artifacts preserved')
 fd=os.open(out,os.O_CREAT|os.O_EXCL|os.O_RDWR,0o600);os.ftruncate(fd,size);chunk=32*1024*1024
 def chunk_read(start):
  end=min(size,start+chunk)-1
  for attempt in range(6):
   try:
    with requests.get(url,headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=(20,90)) as r:
     r.raise_for_status();require(r.status_code==206 and r.headers['Content-Range']==f'bytes {start}-{end}/{size}','Source range identity differs')
     position=start
     for b in r.iter_content(1024*1024):
      require(position+len(b)<=end+1,'Source exceeded requested range');require(os.pwrite(fd,b,position)==len(b),'Short staging write');position+=len(b)
     require(position==end+1,'Truncated source range')
    return
   except Exception:
    if attempt==5:raise
    time.sleep(2**attempt)
 try:
  with concurrent.futures.ThreadPoolExecutor(4) as pool:list(pool.map(chunk_read,range(0,size,chunk)))
  os.fsync(fd)
 finally:os.close(fd)
 actual=sha(out);require(not spec.get('sha256') or actual==spec['sha256'],'Source SHA256 differs')
 if spec.get('md5'):
  h=hashlib.md5()
  with out.open('rb') as f:
   while b:=f.read(4*1024*1024):h.update(b)
  require(h.hexdigest()==spec['md5'],'Source MD5 differs')
 require(spec.get('sha256') or spec.get('md5'),'Expected published source digest required')
 write(out.with_suffix('.receipt.json'),{'source':source_id,'url':spec['url'],'downloadURL':url,'bytes':size,'sha256':actual,'verified':True})
 return str(out)

def sources(cohort,task,output):
 """A separate, relocatable invocation map; never edits a sealed corpus receipt."""
 root=Path(cohort).resolve();require(task in TASKS,'Unknown development task')
 a=read(root/'assay.json')
 manifest={'schemaVersion':1,'format':'numivivo-cell-response-cli-sources/v1','sources':[{'corpus':str(root/s['directory']/('corpus-'+task)),'store':str(root/s['directory']/'store')} for s in a['sources']]}
 write(output,manifest);return output

def metadata(manifest,output,source_id=None,source=None):
 """Read only the source observation/feature metadata, never the expression matrix."""
 import fsspec
 m=read(manifest);out=Path(output);out.mkdir(parents=True,exist_ok=True)
 for spec in m['sources']:
  if source_id and spec['id']!=source_id:continue
  path=out/(spec['id']+'.json');require(not path.exists(),'Metadata already exists: '+str(path))
  origin=Path(source).open('rb') if source else fsspec.open(spec.get('downloadURL',spec['url']),block_size=1048576).open()
  with origin as stream,h5py.File(stream,'r') as h:
   x=h[spec['matrixPath']]
   value={'id':spec['id'],'url':spec['url'],'obs':{k:column(h['obs'],k) for k in ['cell_barcode','gem_group','gene','sgID_AB']},'var':{k:column(h['var'],k) for k in ['gene_id','gene_name']},'matrix':{'shape':list(x.shape) if isinstance(x,h5py.Dataset) else list(x.attrs['shape'])}}
  write(path,value)
 return str(out)

def freeze(manifest,metadata,output):
 m=read(manifest);out=Path(output);require(not out.exists(),'Registration already exists')
 require((m['support']['targetLimit'] is None or 1<=m['support']['targetLimit']<=256) and m['support']['minimumCells']>=32,'Admission support bounds')
 sources=[];counts={};features=set()
 for s in m['sources']:
  p=Path(metadata)/(s['id']+'.json');data=read(p);require(data['id']==s['id'],'Source identity differs')
  c=collections.Counter(data['obs']['gene']);counts[s['id']]=c;features.update(data['var']['gene_id'])
  require(len(set(data['obs']['cell_barcode']))==len(data['obs']['cell_barcode']),'Source barcodes not unique')
  require(c[s['controlLabel']]>=m['support']['minimumControls'],'Insufficient declared controls')
  require(len(set(data['var']['gene_id']))==len(data['var']['gene_id']),'Ambiguous feature axis')
  sources.append({**s,'metadataSHA256':sha(p),'sourceCells':len(data['obs']['gene']),'sourceFeatures':len(data['var']['gene_id'])})
 targets=sorted(set().union(*(set(c) for c in counts.values()))-{s['controlLabel'] for s in sources})
 eligible=[t for t in targets if sum(c[t]>=m['support']['minimumCells'] for c in counts.values())>=m['support']['minimumContexts']]
 selected=sorted(sorted(eligible,key=lambda t:digest(m['seed']+'|target|'+t))[:m['support']['targetLimit']])
 held=sorted(selected,key=lambda t:digest(m['seed']+'|held-target|'+t))[:max(1,len(selected)//5)]
 registration={**m,'sources':sources,'targets':selected,'heldTargets':held,'featureAxis':sorted(features),'metadataOnlySelection':True,'biologicalPromotion':False,
  'exclusions':[{'target':t,'reason':'below metadata support' if t not in eligible else 'outside registered deterministic engineering limit','counts':{sid:c[t] for sid,c in counts.items()}} for t in targets if t not in selected],
  'splitScope':'Technical development transfer only. GEM groups are captures; biological unit and donor independence unresolved.',
  'sourcePanelLimitation':'Deposits prefilter their gene panels using source-wide expression. No new expression-based feature selection is performed; absent source features are unmeasured, not zeros.',
  'modelSelection':'No fitting or checkpoint selection in admission. Training, target representations, and comparison budgets require a separate registration.'}
 write(out,registration);return registration

def partitions(r,s,target,block,control=False):
 block_partition='validation' if int(digest(r['seed']+'|block|'+s['id']+'|'+str(block))[:8],16)%5==0 else 'training'
 return {'known-target':block_partition,'held-target':'validation' if not control and target in r['heldTargets'] else 'training','held-context':'validation' if s['id'] in r['heldContexts'] else 'training'}

def row_selection(r,s,m):
 obs=m['obs'];targets=set(r['targets']);rows=[]
 for target in sorted(targets|{s['controlLabel']}):
  ids=[i for i,t in enumerate(obs['gene']) if t==target]
  if target!=s['controlLabel'] and len(ids)<r['support']['minimumCells']:continue
  limit=r['support']['controlCellLimit'] if target==s['controlLabel'] else r['support']['cellsPerTarget']
  chosen=sorted(ids,key=lambda i:digest(r['seed']+'|cell|'+s['id']+'|'+obs['cell_barcode'][i]))[:limit] if limit else ids
  rows.extend(chosen)
 result=[]
 for i in sorted(rows):
  target=obs['gene'][i];control=target==s['controlLabel'];block=str(obs['gem_group'][i]);parts=partitions(r,s,target,block,control)
  result.append({'sourceRow':i,'id':s['id']+':'+obs['cell_barcode'][i],'barcode':obs['cell_barcode'][i],'source':s['id'],'context':s['context'],'experimentalBlock':s['id']+':gem-'+block,'biologicalUnit':None,'biologicalUnitStatus':'unresolved; capture is not a biological replicate','guideAssignment':obs['sgID_AB'][i] or None,'target':None if control else target,'control':control,'interventionType':'CRISPRi','measurementAxis':s['id'],'partitions':parts,'observationAccess':'development-open' if control or all(p=='training' for p in parts.values()) else 'reserved-development'})
 return result

def native(binary,*args,output):
 p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True)
 require(p.returncode==0,p.stderr[-6000:]);Path(output).write_text(p.stdout);return json.loads(p.stdout)

def admit(registration,metadata,source_id,source,output,binary):
 r=read(registration);s=next(x for x in r['sources'] if x['id']==source_id);meta=read(Path(metadata)/(source_id+'.json'));require(sha(Path(metadata)/(source_id+'.json'))==s['metadataSHA256'],'Metadata changed after registration')
 source=Path(source);out=Path(output)/source_id;out.mkdir(parents=True,exist_ok=False)
 # Verify whole source before any count processing. Metadata and selected bytes
 # are subsequently compared against the same source, not a different mirror.
 actual=sha(source);require(not s.get('sha256') or actual==s['sha256'],'Source SHA256 mismatch')
 if s.get('md5'):
  h=hashlib.md5()
  with source.open('rb') as f:
   while b:=f.read(4*1024*1024):h.update(b)
  require(h.hexdigest()==s['md5'],'Source MD5 mismatch')
 rows=row_selection(r,s,meta);require(rows,'No supported rows')
 write(out/'rows.json',rows);write(out/'source.json',{**s,'verifiedSHA256':actual,'registrationSHA256':sha(registration),'projectionRowsSHA256':sha(out/'rows.json')})
 ids=[x['sourceRow'] for x in rows];features=meta['var']['gene_id'];names=meta['var']['gene_name']
 # Write CSR incrementally: one source row and at most 256 sparse rows in memory.
 import anndata,pandas as pd
 samples={};sample_ids=[]
 for cell in rows:
  target=cell['target'] or s['controlLabel'];part=cell['partitions']['known-target'];sid=s['id']+'|'+target+'|'+part
  sample_ids.append(sid);samples[sid]={'id':sid,'biologicalReplicateID':'unresolved:'+s['id'],'condition':target,'batchID':s['context'],'organism':'Homo sapiens'}
 a=anndata.AnnData(sparse.csr_matrix((len(rows),len(features)),dtype=np.uint64),obs=pd.DataFrame({'sample':sample_ids,'barcode':[c['id'] for c in rows]},index=[c['id'] for c in rows]),var=pd.DataFrame({'feature_id':features,'feature_name':names},index=features))
 a.write_h5ad(out/'projection.h5ad',compression='gzip');del a
 row_checks=[];indptr=[0];buffer_values=[];buffer_indices=[];written=0
 with h5py.File(source,'r') as h,h5py.File(out/'projection.h5ad','a') as projection:
  for k in ('cell_barcode','gene','gem_group','sgID_AB'):require(column(h['obs'],k)==meta['obs'][k],'Metadata/source mismatch: '+k)
  require(column(h['var'],'gene_id')==features,'Feature axis changed')
  x=h[s['matrixPath']];dest=projection['X']
  for key,dtype in [('data','<u8'),('indices','<i4')]:
   del dest[key];dest.create_dataset(key,shape=(0,),maxshape=(None,),dtype=dtype,chunks=(65536,),compression='gzip')
  for n,i in enumerate(ids):
   if isinstance(x,h5py.Dataset):v=np.asarray(x[i]);ix=np.flatnonzero(v);v=v[ix]
   else:
    require(x.attrs.get('encoding-type','csr_matrix') in ('csr_matrix',b'csr_matrix'),'CSR or row-major dense source required')
    lo,hi=map(int,x['indptr'][i:i+2]);ix=x['indices'][lo:hi];v=x['data'][lo:hi]
   require(np.all(np.isfinite(v)) and np.all(v>=0) and np.all(v==np.floor(v)) and np.all(v<2**64),'Raw count source contains unsupported values')
   require(v.sum()>0,'Empty selected source cell; retain failure for preregistration revision')
   v=v.astype(np.uint64);ix=ix.astype(np.int32);buffer_values.append(v);buffer_indices.append(ix);indptr.append(indptr[-1]+len(v))
   if n%max(1,len(ids)//16)==0:row_checks.append({'row':n,'sourceRow':i,'nonzeros':len(v),'sha256':hashlib.sha256(ix.astype('<u4').tobytes()+v.astype('<u8').tobytes()).hexdigest()})
   if len(buffer_values)==256 or n==len(ids)-1:
    stop=indptr[-1]
    for key,b in [('data',buffer_values),('indices',buffer_indices)]:dest[key].resize((stop,));dest[key][written:stop]=np.concatenate(b)
    written=stop;buffer_values=[];buffer_indices=[]
  del dest['indptr'];dest.create_dataset('indptr',data=np.asarray(indptr,dtype='<i8'))
 mapping={'schemaVersion':1,'id':s['id'],'evidence':'measured','sourceDescription':s['url']+'; exact registered engineering rows; see source.json','countUnit':'umiCount','matrixPath':'X','samples':list(samples.values()),'sampleColumn':'sample','barcodeColumn':'barcode','featureIDColumn':'feature_id','featureNameColumn':'feature_name','mitochondrialFeatureIDs':[]}
 write(out/'mapping.json',mapping)
 native(binary,'singlecell-h5ad-store',out/'projection.h5ad','--plan',out/'mapping.json','--output',out/'store',output=out/'native-store.json')
 # Native records are source-major tuples. Compare sampled rows exactly, and
 # bind all bytes through the native receipt; no float-count reconstruction.
 records=np.memmap(out/'store/counts.bin',dtype=[('row','<u4'),('feature','<u4'),('count','<u8')],mode='r');offsets=np.r_[0,np.cumsum(np.bincount(records['row'],minlength=len(rows)))]
 require(np.all(records['row'][1:]>=records['row'][:-1]),'Native count store row order is unsupported')
 for c in row_checks:
  v=records[offsets[c['row']]:offsets[c['row']+1]];require(hashlib.sha256(v['feature'].tobytes()+v['count'].tobytes()).hexdigest()==c['sha256'],'Native/source count reconstruction differs')
 np.save(out/'row-offsets.npy',offsets);write(out/'row-reconstruction.json',{'exact':True,'samples':row_checks,'nativeReceiptSHA256':sha(out/'store/receipt.json')})
 # An explicit structural descriptor enables corpus admission only. It is NOT a
 # biological representation and cannot qualify a learned target comparison.
 descriptor={'schemaVersion':1,'format':'numivivo-cell-response-descriptor-source/v1','sourceDescription':'Admission-only constant descriptor; replace under a training preregistration before fitting. No biological target information.','targets':[{'id':t,'descriptors':[0]} for t in r['targets']],'sourceArtifacts':[fingerprint(sha(registration))]}
 (out/'descriptors.json').write_text(json.dumps(descriptor,sort_keys=True,separators=(',',':'),ensure_ascii=False))
 descriptor_hash=sha(out/'descriptors.json')
 for task in TASKS:
  assignments=[]
  # Source samples retain technical-block partition. For held-target validation,
  # dedicate the same registered control partition; controls cannot be reused
  # as distinct cells or distinct biological units.
  for sid,sample in samples.items():
   cells=[c for c,ss in zip(rows,sample_ids) if ss==sid];first=cells[0];iscontrol=first['control'];part=first['partitions'][task]
   if task=='held-target' and iscontrol:part=first['partitions']['known-target']
   assignments.append({'sampleID':sid,'sourceSample':sample,'role':'control' if iscontrol else 'perturbed','targetID':None if iscontrol else first['target'],'guideID':sample['condition'],'modality':'CRISPRi','studyID':s['study'],'contextID':s['context'],'pairID':s['id'],'partition':part})
  for x in assignments:
   require(any(y['role']=='control' and y['partition']==x['partition'] for y in assignments),'No partition-matched controls; revise source support rather than duplicate controls')
  plan={'schemaVersion':3,'id':s['id']+'-'+task,'sourceDescription':'Technical development '+task+'; exact guides/blocks in rows.json. Legacy guideID stores target label only.','featureAxis':{'featureIDs':r['featureAxis'],'normalizationTarget':10000,'missingValue':-1},'targets':descriptor['targets'],'descriptorSource':fingerprint(descriptor_hash),'assignments':assignments,'heldOutTargetIDs':r['heldTargets'] if task=='held-target' else [],'heldOutContextIDs':[x['context'] for x in r['sources'] if x['id'] in r['heldContexts']] if task=='held-context' else []}
  write(out/(task+'.json'),plan)
  native(binary,'cell-response-prepare',out/'store','--plan',out/(task+'.json'),'--descriptor-source',out/'descriptors.json','--output',out/('corpus-'+task),output=out/(task+'-receipt.json'))
 write(out/'admission.json',{'source':s['id'],'cells':len(rows),'targets':len({c['target'] for c in rows if not c['control']}),'controls':sum(c['control'] for c in rows),'features':len(features),'rawCounts':True,'sparseEntries':len(records),'biologicalUnits':None,'guideIdentity':'source sgID_AB retained per row; no inference','normalization':'Read time log1p(count / retained-source-panel cell total * 10000); raw counts retained','sourcePanelLimit':r['sourcePanelLimitation']})
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 sm=sp.add_parser('sources');sm.add_argument('--cohort',required=True);sm.add_argument('--task',choices=TASKS,required=True);sm.add_argument('--output',required=True)
 ft=sp.add_parser('fetch');ft.add_argument('--manifest',required=True);ft.add_argument('--source-id',required=True);ft.add_argument('--output',required=True)
 mt=sp.add_parser('metadata');mt.add_argument('--manifest',required=True);mt.add_argument('--output',required=True);mt.add_argument('--source-id');mt.add_argument('--source')
 f=sp.add_parser('freeze');f.add_argument('--manifest',required=True);f.add_argument('--metadata',required=True);f.add_argument('--output',required=True)
 a=sp.add_parser('admit');a.add_argument('--registration',required=True);a.add_argument('--metadata',required=True);a.add_argument('--source-id',required=True);a.add_argument('--source',required=True);a.add_argument('--output',required=True);a.add_argument('--binary',required=True)
 args=vars(p.parse_args());cmd=args.pop('cmd');print({'freeze':freeze,'admit':admit,'metadata':metadata,'fetch':fetch,'sources':sources}[cmd](**args))
