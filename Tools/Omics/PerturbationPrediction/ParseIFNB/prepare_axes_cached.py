from remote_io import *
import numpy as np, concurrent.futures, sys, os
ROSTER=json.loads((ROOT/'sources/roster.json').read_text());CODES=np.load(ROOT/'sources/roster-codes.npz');selected=CODES['selected'];cats=ROSTER['categories']
starts=np.r_[0,np.flatnonzero(np.diff(selected)!=1)+1];ends=np.r_[starts[1:],len(selected)]
# The original HDF5 metadata lies near its end. This is an ephemeral RAM cache,
# never a substitute for ETag/range checks and never written as a source snapshot.
BASE=225900000000; STEP=16000000
cache={};ranges=[]
def fetch(start):
 stop=min(SIZE,start+STEP)-1;req=urllib.request.Request(URL,headers={'Range':f'bytes={start}-{stop}','If-Match':ETAG,'Accept-Encoding':'identity'});before=time.time()
 with urllib.request.urlopen(req,timeout=120) as response:
  assert response.status==206 and response.headers['Content-Range']==f'bytes {start}-{stop}/{SIZE}' and response.headers['ETag']==ETAG
  data=response.read(STEP+1);assert len(data)==stop-start+1
 return start,data,dict(start=start,stop=stop,SHA256=hashlib.sha256(data).hexdigest(),seconds=time.time()-before)
class Cached(Remote):
 def read(self,n=-1):
  if self.position>=BASE and n>=0:
   end=min(SIZE,self.position+n);parts=[]
   while self.position<end:
    start=BASE+(self.position-BASE)//STEP*STEP;data=cache[start];offset=self.position-start;take=min(end-self.position,len(data)-offset);parts.append(data[offset:offset+take]);self.position+=take
   return b''.join(parts)
  return super().read(n)
if __name__=='__main__':
 out=ROOT/'axes';out.mkdir(exist_ok=True)
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
   for start,data,receipt in executor.map(fetch,range(BASE,SIZE,STEP)):
    cache[start]=data;ranges.append(receipt);print('RAM-CACHE',len(cache),sum(len(x) for x in cache.values()),flush=True)
  remote=Cached()
  with h5py.File(remote,'r',rdcc_nbytes=0) as h:
   for i,(lo,hi) in enumerate(zip(starts,ends)):
    if (out/f'run-{i:04d}.npz').exists():continue
    start,end=int(selected[lo]),int(selected[hi-1])+1
    ptr=h['X/indptr'][start:end+1];barcodes=h['obs/_index'].asstr()[start:end]
    assert len(ptr)==end-start+1 and (np.diff(ptr)>=0).all() and 0<=ptr[0]<=ptr[-1]<=18830591942
    totals=h['obs/tscp_count'][start:end];genes=h['obs/gene_count'][start:end]
    np.savez_compressed(out/f'run-{i:04d}.npz',lo=int(lo),hi=int(hi),sourceRows=selected[lo:hi],indptr=ptr,barcodes=barcodes.astype(str),tscp_count=totals,gene_count=genes)
    if i%100==0:print('AXIS',i,int(lo),int(hi),flush=True)
  ranges+=remote.requests
  print('AXES COMPLETE',len(starts),flush=True)
 finally:(out/'cached-metadata-ranges.json').write_text(json.dumps(dict(url=URL,etag=ETAG,size=SIZE,requests=ranges,bytesRead=sum(r['stop']-r['start']+1 for r in ranges)),indent=2)+'\n')
