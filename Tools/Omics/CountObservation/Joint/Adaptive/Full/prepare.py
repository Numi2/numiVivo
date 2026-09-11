"""Freeze full-gene sparse training cells, one donor/condition CSR group at a time.

No dense cells by genes array. The CSC cache retains every selected count and
original full RNA denominator. Counts and row/feature identities are checked
against the original hashed H5AD and the full published calibration.
"""
from pathlib import Path
import collections,gzip,hashlib,json,os,sys,time,resource
import h5py,numpy as np
from scipy import sparse
root=Path(sys.argv[1]);origin=sys.argv[2];parent=Path('/Users/home/numivivo-count-calibration-20260911')
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(8<<20):h.update(b)
 return h.hexdigest()
meta=read(parent/'inputs'/(origin+'.json'));source=Path(meta['sourcePath']);report=read(parent/(origin+'-native.json.gz'))
assert report['metadataSHA256']==sha(parent/'inputs'/(origin+'.json'))
state={'status':'running','origin':origin,'pid':os.getpid(),'startedUnix':time.time(),'sourceCells':0,'nonzeros':0,'groups':[]}
def save(): (root/(origin+'-prepare-state.json')).write_text(json.dumps(state,indent=2)+'\n')
save();assert sha(source)==meta['sourceSHA256'];st=source.stat();signature=(st.st_ino,st.st_size,st.st_mtime_ns)
output=root/(origin+'-cells.h5');temp=root/(origin+'-cells.partial.h5');assert not output.exists() and not temp.exists()
try:
 with h5py.File(source,'r') as h,h5py.File(temp,'x') as out:
  assert h['X'].attrs['shape'].tolist()==meta['sourceShape'];indptr=h['X/indptr'][:];gcount=len(meta['featureIDs'])
  out.attrs['schema']='numivivo-full-gene-training-csc-v1';out.attrs['sourceSHA256']=meta['sourceSHA256'];out.attrs['metadataSHA256']=sha(parent/'inputs'/(origin+'.json'))
  out.create_dataset('featureIDs',data=np.array(meta['featureIDs'],dtype=h5py.string_dtype('utf-8')))
  for gi,group in enumerate(meta['groups']):
   rows=np.array(sorted(group['sourceRows']),np.int64);nnz_per=indptr[rows+1]-indptr[rows];assert np.all(nnz_per>0)
   ptr=np.concatenate(([0],np.cumsum(nnz_per))).astype(np.int64);nnz=int(ptr[-1]);ix=np.empty(nnz,np.int32);yy=np.empty(nnz,np.uint32);start=0
   while start<len(rows):
    end=start+1
    while end<len(rows) and end-start<512 and rows[end]==rows[end-1]+1:end+=1
    first,last=int(indptr[rows[start]]),int(indptr[rows[end-1]+1]);a,b=int(ptr[start]),int(ptr[end]);assert b-a==last-first
    values=h['X/data'][first:last];assert np.all(np.isfinite(values)) and np.all(values>0) and np.all(values==np.floor(values)) and np.max(values)<=1_000_000_000
    ix[a:b]=h['X/indices'][first:last];yy[a:b]=values;start=end
   matrix=sparse.csr_matrix((yy,ix,ptr),shape=(len(rows),gcount));assert matrix.has_canonical_format and matrix.has_sorted_indices
   libraries=np.add.reduceat(yy.astype(np.uint64),ptr[:-1]);assert np.all(libraries>0) and np.max(libraries)<=1_000_000_000
   counts=np.bincount(ix,weights=yy,minlength=gcount).astype(np.uint64);previous=report['moments'][gi]
   assert counts.tolist()==previous['counts'] and int(libraries.sum())==previous['libraryCounts'] and len(rows)==previous['cells']
   csc=matrix.tocsc();roundtrip=csc.tocsr();assert np.array_equal(roundtrip.indptr,ptr) and np.array_equal(roundtrip.indices,ix) and np.array_equal(roundtrip.data,yy);del roundtrip
   groupout=out.create_group(str(gi));groupout.attrs['donorID']=group['donorID'];groupout.attrs['conditionID']=group['conditionID'];groupout.attrs['shape']=[len(rows),gcount]
   depth,depth_index,depth_cells=np.unique(libraries,return_inverse=True,return_counts=True)
   arrays={'sourceRows':rows,'libraryCounts':libraries.astype(np.uint32),'depthLibraries':depth.astype(np.uint32),'depthIndex':depth_index.astype(np.uint32),'depthCells':depth_cells.astype(np.uint32),'indptr':csc.indptr.astype(np.uint64),
    'indices':csc.indices.astype(np.uint16 if len(rows)<=65535 else np.uint32),'data':csc.data.astype(np.uint16 if csc.data.max()<=65535 else np.uint32)}
   hashes={}
   for name,array in arrays.items():
    ds=groupout.create_dataset(name,data=array,compression='gzip',compression_opts=4,shuffle=True,chunks=True)
    # Exact readback, including original integer values after compact storage.
    assert np.array_equal(ds[:],array);hashes[name]=hashlib.sha256(array.tobytes()).hexdigest()
   out.flush();state['sourceCells']+=len(rows);state['nonzeros']+=nnz;state['groups'].append({'index':gi,'donorID':group['donorID'],'conditionID':group['conditionID'],'cells':len(rows),'nonzeros':nnz,'depthBins':len(depth),'arraySHA256':hashes});state['cacheBytesSoFar']=temp.stat().st_size;state['peakRSSBytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;save();print(json.dumps({k:v for k,v in state.items() if k!='groups'}),flush=True)
   del matrix,csc,arrays,ix,yy,ptr,libraries,counts
 assert signature==(source.stat().st_ino,source.stat().st_size,source.stat().st_mtime_ns)
 assert state['sourceCells']==report['cells'] and state['nonzeros']==report['nonzeros']
 temp.replace(output)
 models={m['conditionID']:{f['featureID']:f for f in m['features']} for m in report['models']};genes=[]
 for j,feature in enumerate(meta['featureIDs']):
  c=models['control'][feature];t=models['IFNB'][feature];genes.append({'featureIndex':j,'featureID':feature,'controlCellDispersion':c.get('cellDispersion'),'treatedCellDispersion':t.get('cellDispersion'),'controlCalibrationStatus':c['status'],'treatedCalibrationStatus':t['status']})
 manifest={'schemaVersion':1,'origin':origin,'source':meta,'sourceCalibrationSHA256':sha(parent/(origin+'-native.json.gz')),'cacheSHA256':sha(output),'genes':genes,'groups':state['groups'],'qualification':'All original admitted cells and full RNA counts, full source gene axis. Fixed previously estimated gene/condition dispersion. No query-treated outcomes or donor-omission refitting.'}
 raw=json.dumps(manifest,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
 with (root/(origin+'-manifest.json.gz')).open('xb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(raw)
 state.update(status='completed',genes=len(genes),availablePairedGenes=sum(g['controlCellDispersion'] is not None and g['treatedCellDispersion'] is not None for g in genes),cacheSHA256=manifest['cacheSHA256'],cacheBytes=output.stat().st_size,manifestSHA256=hashlib.sha256(raw).hexdigest())
except BaseException as e:state.update(status='failed',error=repr(e));raise
finally:state.update(finishedUnix=time.time(),seconds=time.time()-state['startedUnix'],peakRSSBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);save()
