"""Stream original sparse cells to the physical Mac native owner; never retain a dense cell matrix."""
from pathlib import Path
import gzip,hashlib,json,os,struct,subprocess,sys,threading,time
import numpy as np
import h5py
root=Path(sys.argv[1]);origin=sys.argv[2];metadata=root/'inputs'/(origin+'.json');rawmeta=metadata.read_bytes();meta=json.loads(rawmeta)
remote='/Users/n/numivivo-count-calibration-20260911';statefile=root/(origin+'-pipeline.json')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
assert sha(metadata)==json.loads((root/'inputs/freeze.json').read_bytes())['files'][metadata.name]
source=Path(meta['sourcePath']);assert sha(source)==meta['sourceSHA256']
source_stat=(source.stat().st_size,source.stat().st_mtime_ns,source.stat().st_ino)
output=root/(origin+'-native.json.gz');assert not output.exists()
subprocess.run(['ssh','-4','-o','ConnectTimeout=10','macmini','mkdir -p '+remote+'/inputs'],check=True)
subprocess.run(['scp','-q',str(metadata),'macmini:'+remote+'/inputs/'+metadata.name],check=True)
state={'status':'running','pid':os.getpid(),'origin':origin,'startedUnix':time.time(),'sourceSHA256':meta['sourceSHA256'],'metadataSHA256':sha(metadata),'cells':0,'nonzeros':0,'sourceGroupsChecked':0}
def save():statefile.write_text(json.dumps(state,indent=2)+'\n')
save();reference=[];stream=hashlib.sha256();binding=hashlib.sha256(rawmeta)
err=(root/(origin+'-native.stderr')).open('wb')
command=['ssh','-C','-4','-o','ConnectTimeout=10','-o','ServerAliveInterval=10','-o','ServerAliveCountMax=1','macmini','/usr/bin/time -l '+remote+'/build/calibrate-cells '+remote+'/inputs/'+metadata.name]
proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err);state['sshPID']=proc.pid;state['command']=command;save()
reader_errors=[]
def collect():
 try:
  with output.open('wb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:
   while block:=proc.stdout.read(1024*1024):g.write(block)
 except BaseException as e:reader_errors.append(repr(e))
reader=threading.Thread(target=collect);reader.start()
try:
 with h5py.File(source,'r') as h:
  ptr=h['X/indptr'][:];m=len(meta['featureIDs']);fullshape=meta['sourceShape'];assert len(ptr)==fullshape[0]+1
  assert h['X'].attrs['shape'].tolist()==fullshape
  assert np.all(ptr[1:]>=ptr[:-1]) and ptr[0]==0 and ptr[-1]==len(h['X/data'])==len(h['X/indices'])
  original=np.load('/Users/home/numivivo-cross-study-ifnb-20260911/inputs/'+origin+'-source-counts.npz')
  assert original['featureIDs'].tolist()==meta['featureIDs']
  for gi,group in enumerate(meta['groups']):
   rows=sorted(group['sourceRows']);n=len(rows);s1=np.zeros(m);s2=np.zeros(m);shot=np.zeros(m);counts=np.zeros(m,np.int64);positive=np.zeros(m,np.int64);library=0;library2=0.0
   # Run-local blocks never include unselected rows or create a dense cell matrix.
   start=0
   while start<n:
    end=start+1
    while end<n and end-start<512 and rows[end]==rows[end-1]+1:end+=1
    row0=rows[start];row1=rows[end-1]+1;first=int(ptr[row0]);last=int(ptr[row1]);offsets=ptr[row0:row1+1]-first
    indices=h['X/indices'][first:last].astype(np.int64);data=h['X/data'][first:last].astype(np.int64)
    assert np.all(indices>=0) and np.all(indices<m) and np.all(data>0) and np.all(data<=np.iinfo(np.uint32).max)
    for k,row in enumerate(range(row0,row1)):
     lo,hi=map(int,offsets[k:k+2]);ii=indices[lo:hi];yy=data[lo:hi]
     assert len(ii)>0 and np.all(ii[1:]>ii[:-1])
     depth=int(yy.sum());assert 0<depth<=1_000_000_000
     z=yy*(1e6/depth);counts[ii]+=yy;positive[ii]+=1;s1[ii]+=z;s2[ii]+=z*z;shot[ii]+=yy*(1e6/depth)**2;library+=depth;library2+=float(depth)**2
     payload=np.empty((len(ii),2),dtype='<u4');payload[:,0]=ii;payload[:,1]=yy
     packet=struct.pack('<III',gi,row,len(ii))+payload.tobytes()
     proc.stdin.write(packet);stream.update(packet);binding.update(packet)
     state['cells']+=1;state['nonzeros']+=len(ii)
    start=end
   np.testing.assert_array_equal(counts,original['counts'][gi]);state['sourceGroupsChecked']+=1
   mean=s1/n;variance=np.maximum(0,(s2-s1*s1/n)/(n-1));pair=np.maximum(0,(s1*s1-s2)/(n*(n-1)));pair[positive<2]=0
   cov=(1e6*counts-library*mean)/(n-1);depthvar=max(0,(library2-library*library/n)/(n-1))
   corr=np.divide(cov,np.sqrt(depthvar*variance),out=np.zeros(m),where=(depthvar*variance)>0)
   reference.append(dict(cells=n,libraryCounts=library,counts=counts.tolist(),positiveCells=positive.tolist(),meanCPM=mean.tolist(),sampleVarianceCPM=variance.tolist(),meanPoissonVarianceCPM=(shot/n).tolist(),distinctCellRateProductCPM2=pair.tolist(),depthRateCorrelation=corr.tolist()))
   save();print(json.dumps({'origin':origin,'groups':gi+1,'cells':state['cells'],'nonzeros':state['nonzeros']}),flush=True)
 proc.stdin.close();state['nativeExit']=proc.wait();reader.join();err.close()
 assert not reader_errors,reader_errors
 assert state['nativeExit']==0
 native=json.loads(gzip.decompress(output.read_bytes()))
 assert native['streamSHA256']==stream.hexdigest() and native['bindingSHA256']==binding.hexdigest()
 assert native['metadataSHA256']==sha(metadata) and native['cells']==state['cells'] and native['nonzeros']==state['nonzeros']
 assert source_stat==(source.stat().st_size,source.stat().st_mtime_ns,source.stat().st_ino),'source changed during scan'
 with (root/(origin+'-reference-moments.json.gz')).open('wb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(json.dumps(reference,sort_keys=True,allow_nan=False).encode())
 state['streamSHA256']=stream.hexdigest();state['bindingSHA256']=binding.hexdigest();state['nativeReportSHA256']=sha(output);state['status']='completed-source-linked-native-calibration-reference-comparison-pending'
except BaseException as e:
 state['status']='failed';state['error']=repr(e)
 if proc.poll() is None:proc.stdin.close();state['nativeExit']=proc.wait()
 reader.join();err.close();raise
finally:state['finishedUnix']=time.time();state['seconds']=state['finishedUnix']-state['startedUnix'];save()
