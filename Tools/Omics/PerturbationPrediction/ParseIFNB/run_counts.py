"""Replay every selected count through the real native owner; never store full X."""
from pathlib import Path
import numpy as np, json, hashlib, http.client, urllib.parse, ssl, threading, time, concurrent.futures, subprocess, sys, resource, os
ROOT=Path(os.environ['NUMIVIVO_PARSE_STUDY']);identity=json.loads((ROOT/'sources/asset-head.json').read_text());URL=identity['url'];ETAG=identity['headers']['ETag'];SIZE=int(identity['headers']['Content-Length']);PARSED=urllib.parse.urlsplit(URL)
RUNS=json.loads((ROOT/'prepared/runs.json').read_text());VERIFY='--verify' in sys.argv;NAME='replay' if VERIFY else 'ingest';OUT=ROOT/NAME;LOCAL=threading.local();FEATURES=40352;MAX_RANGE=16000000

def fetch(start,stop,receipts):
 # Retry only transport failures, with the same immutable If-Match condition.
 errors=[]
 for attempt in range(6):
  if not hasattr(LOCAL,'connection'):LOCAL.connection=http.client.HTTPSConnection(PARSED.netloc,timeout=90,context=ssl.create_default_context())
  before=time.monotonic()
  try:
   LOCAL.connection.request('GET',PARSED.path,headers={'Range':f'bytes={start}-{stop}','If-Match':ETAG,'Accept-Encoding':'identity'})
   response=LOCAL.connection.getresponse()
   if response.status in {500,502,503,504}:
    body=response.read(4096);raise http.client.HTTPException(f'Retryable source HTTP {response.status}; body prefix SHA256 {hashlib.sha256(body).hexdigest()}')
   assert response.status==206 and response.getheader('Content-Range')==f'bytes {start}-{stop}/{SIZE}' and response.getheader('ETag')==ETAG,(response.status,response.getheader('Content-Range'))
   data=response.read(stop-start+2);assert len(data)==stop-start+1
   receipts.append(dict(start=start,stop=stop,SHA256=hashlib.sha256(data).hexdigest(),seconds=time.monotonic()-before,transportErrors=errors));return data
  except (OSError,http.client.HTTPException) as error:
   errors.append(dict(type=type(error).__name__,message=str(error)));LOCAL.connection.close();del LOCAL.connection
   if attempt==5:raise
   time.sleep(min(16,2**attempt))
 raise AssertionError('unreachable')

def read_array(spec,start,end,receipts):
 step=spec['step'];dtype=np.dtype(spec['dtype']);chunks=spec['chunks'];result=np.empty(end-start,dtype=dtype)
 assert chunks and len(chunks)==(end+step-1)//step-start//step
 low=min(x[1] for x in chunks);high=max(x[1]+x[2] for x in chunks)
 assert high-low<=128000000,'source run exceeds bounded 128 MB window'
 data=b''.join(fetch(a,min(high,a+MAX_RANGE)-1,receipts) for a in range(low,high,MAX_RANGE))
 covered=0
 for origin,offset,size in chunks:
  left=max(start,origin);right=min(end,origin+step);assert left<right
  values=np.frombuffer(data,dtype=dtype,count=(right-left),offset=offset-low+(left-origin)*dtype.itemsize)
  result[left-start:right-start]=values;covered+=right-left
 assert covered==end-start
 return result

def one(run):
 i=run['run'];axis=ROOT/f'axes/run-{i:04d}.npz';assert hashlib.sha256(axis.read_bytes()).hexdigest()==run['axisSHA256'];a=np.load(axis);spec=json.loads((ROOT/f'chunks/run-{i:04d}.json').read_text());receipts=[];start,end=run['entryFirst'],run['entryLast']
 values=read_array(spec['data'],start,end,receipts);indices=read_array(spec['indices'],start,end,receipts)
 assert np.isfinite(values).all() and (values>0).all() and (values==np.floor(values)).all() and (values<2**32).all()
 assert (indices>=0).all() and (indices<FEATURES).all()
 counts=values.astype(np.uint64);ptr=a['indptr']-start;local_rows=np.repeat(np.arange(run['lo'],run['hi'],dtype=np.uint64),np.diff(ptr));assert len(local_rows)==len(counts)
 # Canonicalize within each source row, rejecting duplicates, never merging them.
 sorted_rows=0
 for row in range(len(ptr)-1):
  lo,hi=int(ptr[row]),int(ptr[row+1]);index=indices[lo:hi]
  if len(index)>1 and (np.diff(index)<=0).any():
   order=np.argsort(index,kind='stable');indices[lo:hi]=index[order];counts[lo:hi]=counts[lo:hi][order];sorted_rows+=1
   assert (np.diff(indices[lo:hi])>0).all(),'duplicate source column'
 packed=np.empty((len(counts),2),dtype='<u8');packed[:,0]=local_rows|(indices.astype(np.uint64)<<np.uint64(32));packed[:,1]=counts
 # Independent exact integer accumulation, separate from Swift's per-record owner.
 bulk=np.zeros(FEATURES,dtype=np.uint64);np.add.at(bulk,indices,counts)
 totals=np.zeros(len(ptr)-1,dtype=np.uint64);positive=np.flatnonzero(np.diff(ptr)>0)
 for row in positive:totals[row]=counts[int(ptr[row]):int(ptr[row+1])].sum(dtype=np.uint64)
 assert int(counts.max(initial=0))*len(counts)<2**64
 qc=dict(run=i,sourceGeneCountMismatch=int(np.count_nonzero(np.diff(ptr)!=a['gene_count'])),sourceTscpCountMismatch=int(np.count_nonzero(totals!=a['tscp_count'])),sourceMinusMatrixTotal=int(a['tscp_count'].sum())-int(totals.sum()),sourceMinusMatrixGeneCount=int(a['gene_count'].sum())-int(np.diff(ptr).sum()),sortedRows=sorted_rows,nonzeros=len(counts),totalCounts=int(totals.sum()),ranges=receipts)
 if VERIFY:
  old=json.loads((ROOT/f'ingest/run-{i:04d}.json').read_text())
  assert [(x['start'],x['stop'],x['SHA256']) for x in receipts]==[(x['start'],x['stop'],x['SHA256']) for x in old['ranges']]
 return run,packed.tobytes(),bulk,totals,qc

if __name__=='__main__':
 assert __import__('shutil').disk_usage(ROOT).free>=650000000, 'Insufficient capacity for retained bundle and reserve'
 OUT.mkdir(exist_ok=False)
 binary=ROOT/'build/numivivo-omics';command=[str(binary),'singlecell-count-stream-verify',str(ROOT/'native')] if VERIFY else [str(binary),'singlecell-count-stream-pseudobulk','--plan',str(ROOT/'prepared/plan.json'),'--output',str(ROOT/'native')]
 groups=sorted(set(r['sample'] for r in RUNS));matrix=np.zeros((len(groups),FEATURES),dtype=np.uint64);totals=np.zeros(725031,dtype=np.uint64);sha=hashlib.sha256();received=0;summary=dict(status='running',command=command,binarySHA256=hashlib.sha256(binary.read_bytes()).hexdigest(),sourceURL=URL,sourceETag=ETAG,sourceBytes=SIZE,mode=NAME,startUTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
 (OUT/'execution.json').write_text(json.dumps(summary,indent=2)+'\n');before=time.monotonic();log=open(OUT/'native.stdout','wb');errorlog=open(OUT/'native.stderr','wb');process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=log,stderr=errorlog)
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
   # A sliding 16-run window bounds queued matrices and stream records.
   pending={};next_submit=0
   for wanted in range(len(RUNS)):
    while next_submit<min(len(RUNS),wanted+16):pending[next_submit]=pool.submit(one,RUNS[next_submit]);next_submit+=1
    run,data,bulk,rowtotals,qc=pending.pop(wanted).result();i=run['run']
    process.stdin.write(data);sha.update(data);received+=len(data);matrix[groups.index(run['sample'])]+=bulk;totals[run['lo']:run['hi']]=rowtotals
    (OUT/f'run-{i:04d}.json').write_text(json.dumps(qc,separators=(',',':'))+'\n')
    if i%50==0:print('STREAMED',i,run['hi'],received,'seconds',time.monotonic()-before,flush=True)
  process.stdin.close();code=process.wait();log.close();errorlog.close();assert code==0,(code,(OUT/'native.stderr').read_text()[-2000:])
  receipt=json.loads((OUT/'native.stdout').read_text());assert receipt['streamBytes']==received
  # VivoFingerprint serializes the exact 32 digest bytes.
  assert bytes(receipt['stream']['bytes']).hex()==sha.hexdigest(),receipt['stream']
  np.savez_compressed(OUT/'independent-counts.npz',groups=np.array(groups),counts=matrix,cellTotals=totals)
  summary.update(status='passed',seconds=time.monotonic()-before,streamBytes=received,streamSHA256=sha.hexdigest(),cells=len(totals),nonzeros=received//16,totalCounts=int(totals.sum()),nativeExit=code,producerPeakRSS=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,nativePeakRSS=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
 except BaseException as error:
  summary.update(status='failed',errorType=type(error).__name__,error=str(error),streamedBytes=received,seconds=time.monotonic()-before)
  if process.stdin and not process.stdin.closed:process.stdin.close()
  process.wait();raise
 finally:
  (OUT/'execution.json').write_text(json.dumps(summary,indent=2)+'\n');log.close();errorlog.close()
 print(json.dumps(summary,indent=2),flush=True)
