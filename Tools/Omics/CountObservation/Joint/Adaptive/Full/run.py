"""Checkpointed full-gene native fitting; bounded input and output gene shards."""
from pathlib import Path
import collections,fcntl,gzip,hashlib,json,os,shlex,subprocess,sys,threading,time
import h5py,numpy as np
root=Path(sys.argv[1]);origin=sys.argv[2];limit=int(sys.argv[3]) if len(sys.argv)>3 else None
workers=int(os.environ.get('NUMIVIVO_FULL_SHARD_COUNT','1'));worker=int(os.environ.get('NUMIVIVO_FULL_SHARD_INDEX','0'));assert 1<=workers<=4 and 0<=worker<workers
out=Path(os.environ.get('NUMIVIVO_FULL_OUTPUT',str(root/origin)));out.mkdir(exist_ok=True);lock=(out/('controller.lock' if workers==1 else f'controller-{worker}-of-{workers}.lock')).open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
state_path=out/('state.json' if workers==1 else f'state-worker-{worker}.json')
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(4<<20):h.update(b)
 return h.hexdigest()
def encoded(v):return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
manifest=read(root/(origin+'-manifest.json.gz'));prepared=read(root/(origin+'-prepare-state.json'));freeze=read(Path(os.environ.get('NUMIVIVO_FULL_RUNTIME',str(root/'runtime-freeze.json'))));binary=freeze['fullBinarySHA256'];cache=root/(origin+'-cells.h5')
assert prepared['status']=='completed' and sha(cache)==manifest['cacheSHA256'];assert hashlib.sha256(gzip.decompress((root/(origin+'-manifest.json.gz')).read_bytes())).hexdigest()==prepared['manifestSHA256']
plan={'binarySHA256':binary,'manifestSHA256':prepared['manifestSHA256'],'shardGenes':64,'shardWorkers':workers}
with (out/'execution-plan.lock').open('a') as plan_lock:
 fcntl.flock(plan_lock,fcntl.LOCK_EX);plan_path=out/'execution-plan.json'
 if plan_path.exists():assert read(plan_path)==plan
 else:plan_path.write_text(json.dumps(plan,indent=2)+'\n')
state={'status':'running','origin':origin,'pid':os.getpid(),'startedUnix':time.time(),'binarySHA256':binary,'manifestSHA256':prepared['manifestSHA256'],'genesTotal':len(manifest['genes']),'workerIndex':worker,'workerCount':workers,'genesCompleted':0,'shardsCompleted':0,'fitStates':{},'activeShard':None}
def save():
 p=state_path;tmp=state_path.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(p)
save();statuses=collections.Counter();started=0
with h5py.File(cache,'r') as h:
 groups=[];depth_indices=[];pointers=[]
 for gi,g in enumerate(manifest['groups']):
  group=h[str(gi)];assert group.attrs['donorID']==g['donorID'] and group.attrs['conditionID']==g['conditionID']
  groups.append({'donorID':g['donorID'],'conditionID':g['conditionID'],'libraryCounts':group['depthLibraries'][:].tolist(),'cellsPerLibrary':group['depthCells'][:].tolist()});depth_indices.append(group['depthIndex'][:]);pointers.append(group['indptr'][:])
 try:
  for shard_index,first in enumerate(range(0,len(manifest['genes']),64)):
   if shard_index%workers!=worker:continue
   last=min(first+64,len(manifest['genes']));tag=f'{first:05d}-{last:05d}';receipt_path=out/(tag+'-receipt.json');input_path=out/(tag+'-input.jsonl.gz');output_path=out/(tag+'-output.jsonl.gz')
   if receipt_path.exists():
    receipt=read(receipt_path);assert receipt['status']=='completed' and receipt['binarySHA256']==binary and receipt['manifestSHA256']==prepared['manifestSHA256'];assert sha(input_path)==receipt['compressedInputSHA256'] and sha(output_path)==receipt['compressedOutputSHA256']
   else:
    if limit is not None and started>=limit:state['status']='bounded-shard-run-completed';break
    started+=1;state['activeShard']=tag;save()
    header={'origin':origin,'firstFeatureIndex':first,'featureIDs':[g['featureID'] for g in manifest['genes'][first:last]],'groups':groups,'cacheSHA256':manifest['cacheSHA256'],'manifestSHA256':prepared['manifestSHA256']}
    if not input_path.exists():
     partial=out/(tag+'-input.partial.gz');assert not partial.exists();input_hash=hashlib.sha256()
     with partial.open('xb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as gz:
      def emit(v):
       data=encoded(v)+b'\n';input_hash.update(data);gz.write(data)
      emit(header)
      for gene in manifest['genes'][first:last]:
       row=dict(gene);row['groups']=[]
       if gene['controlCellDispersion'] is not None and gene['treatedCellDispersion'] is not None:
        j=gene['featureIndex']
        for gi,g in enumerate(groups):
         group=h[str(gi)];lo,hi=map(int,pointers[gi][j:j+2]);rows=group['indices'][lo:hi];values=group['data'][lo:hi]
         counts=np.bincount(depth_indices[gi][rows],weights=values,minlength=len(g['libraryCounts'])).astype(np.uint64);bins=np.flatnonzero(counts)
         row['groups'].append({'bins':bins.tolist(),'counts':counts[bins].tolist()})
       emit(row)
     partial.replace(input_path)
    else:
     input_hash=hashlib.sha256()
     with gzip.open(input_path,'rb') as f:
      while b:=f.read(1<<20):input_hash.update(b)
    native=freeze.get('nativePath','/Users/n/numivivo-full-joint-counts-20260911/full-counts');verify='import hashlib,pathlib;assert hashlib.sha256(pathlib.Path('+repr(native)+').read_bytes()).hexdigest()=='+repr(binary)
    command='python3 -c '+shlex.quote(verify)+' && /usr/bin/time -l '+shlex.quote(native);errors=[];t0=time.time();partial=out/(tag+'-output.partial.gz');assert not partial.exists() and not output_path.exists()
    with (out/(tag+'-stderr.log')).open('xb') as err:
     p=subprocess.Popen(['ssh','-4','-C','-o','ConnectTimeout=10','-o','ServerAliveInterval=10','-o','ServerAliveCountMax=1','macmini',command],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err);state['sshPID']=p.pid;save()
     def send():
      try:
       ih=hashlib.sha256()
       with gzip.open(input_path,'rb') as f:
        while data:=f.read(1<<20):ih.update(data);p.stdin.write(data)
       assert ih.hexdigest()==input_hash.hexdigest()
      except BaseException as e:errors.append(repr(e))
      finally:p.stdin.close()
     writer=threading.Thread(target=send);writer.start();output_hash=hashlib.sha256();output_bytes=0
     with partial.open('xb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0,compresslevel=3) as gz:
      while data:=p.stdout.read(1<<20):gz.write(data);output_hash.update(data);output_bytes+=len(data)
     code=p.wait();writer.join();assert code==0 and not errors,(code,errors)
    shard_states=collections.Counter();count=0
    with gzip.open(partial,'rt') as f:
     for line in f:
      row=json.loads(line);gene=manifest['genes'][first+count];assert row['origin']==origin and row['featureIndex']==gene['featureIndex'] and row['featureID']==gene['featureID'];shard_states[row['status']]+=1;count+=1
    assert count==last-first;partial.replace(output_path)
    receipt={'status':'completed','first':first,'last':last,'genes':count,'fitStates':dict(shard_states),'binarySHA256':binary,'manifestSHA256':prepared['manifestSHA256'],'expandedInputSHA256':input_hash.hexdigest(),'compressedInputSHA256':sha(input_path),'expandedOutputSHA256':output_hash.hexdigest(),'expandedOutputBytes':output_bytes,'compressedOutputSHA256':sha(output_path),'seconds':time.time()-t0}
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
   statuses.update(receipt['fitStates']);state.update(genesCompleted=state['genesCompleted']+receipt['genes'],lastFeatureIndex=last,shardsCompleted=state['shardsCompleted']+1,fitStates=dict(statuses),activeShard=None,updatedUnix=time.time());save();print(json.dumps(state),flush=True)
  else:state['status']='completed-all-source-genes' if workers==1 else 'completed-assigned-source-shards'
 except BaseException as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(finishedUnix=time.time(),seconds=time.time()-state['startedUnix']);save()
