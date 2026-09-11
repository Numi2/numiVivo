"""Cache actual v1 B-tree nodes; HDF5 itself enumerates every chunk afterward.
HDF Group FMT4, section III.A.1: 1D chunk keys include the element-size axis.
The two roots below were observed in HDF5 reads (tree.log), not inferred from
payload positions. This narrowly qualified source helper is not a generic parser.
"""
from remote_io import *
import numpy as np, concurrent.futures, threading, struct
LOCAL=threading.local();CACHE={};RECEIPTS=[];NODE_SIZE=2096;ROOTS={'data':(6848,2),'indices':(75327800048,3)}
def fetch(task):
 address,level=task
 if not hasattr(LOCAL,'remote'):LOCAL.remote=Remote()
 before=time.time();data=LOCAL.remote.fetch(address,address+NODE_SIZE-1,90)
 assert data[:4]==b'TREE' and data[4]==1 and data[5]==level and 0<=level<=3
 count=int.from_bytes(data[6:8],'little');assert 1<=count<=64
 children=[] if level==0 else [(int.from_bytes(data[48+32*i:56+32*i],'little'),level-1) for i in range(count)]
 assert all(0<a<SIZE-NODE_SIZE for a,l in children)
 return address,data,children,dict(start=address,stop=address+NODE_SIZE-1,SHA256=hashlib.sha256(data).hexdigest(),seconds=time.time()-before)
class Cached(Remote):
 def read(self,n=-1):
  data=CACHE.get(self.position)
  if data is not None and n<=len(data):self.position+=n;return data[:n]
  return super().read(n)
if __name__=='__main__':
 frontier=list(ROOTS.values());seen={x[0] for x in frontier};out=ROOT/'sources';before=time.time()
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
   while frontier:
    next_level=[]
    for address,data,children,receipt in pool.map(fetch,frontier):
     CACHE[address]=data;RECEIPTS.append(receipt)
     for child in children:
      assert child[0] not in seen,'cyclic/shared chunk-index node';seen.add(child[0]);next_level.append(child)
     assert len(seen)<10000
     if len(CACHE)%100==0:print('INDEX-NODES',len(CACHE),time.time()-before,flush=True)
    frontier=next_level
  addresses=sorted(CACHE);np.savez_compressed(out/'chunk-index-nodes.npz',addresses=np.array(addresses,dtype=np.uint64),bytes=np.stack([np.frombuffer(CACHE[a],dtype=np.uint8) for a in addresses]))
  f=Cached();maps={}
  with h5py.File(f,'r',rdcc_nbytes=0) as h:
   for key in ['data','indices']:
    d=h['X/'+key];step=d.chunks[0];items=[]
    def visit(info):
     assert info.filter_mask==0 and info.size==step*d.dtype.itemsize
     items.append([int(info.chunk_offset[0]),int(info.byte_offset),int(info.size)])
    d.id.chunk_iter(visit);values=np.array(items,dtype=np.int64);assert np.array_equal(values[:,0],np.arange(len(values))*step)
    assert len(values)==(18830591942+step-1)//step
    maps[key]=dict(step=step,dtype=d.dtype.str,values=values);print('HDF5-ENUMERATED',key,len(values),flush=True)
   np.savez_compressed(out/'all-chunk-map.npz',**{k:v['values'] for k,v in maps.items()})
  chunkdir=ROOT/'chunks';chunkdir.mkdir(exist_ok=True);checked=0
  for i in range(3456):
   a=np.load(ROOT/f'axes/run-{i:04d}.npz');start,end=int(a['indptr'][0]),int(a['indptr'][-1]);result={}
   for key,spec in maps.items():
    step=spec['step'];result[key]=dict(step=step,dtype=spec['dtype'],chunks=spec['values'][start//step:(end+step-1)//step].tolist())
   path=chunkdir/f'run-{i:04d}.json'
   if path.exists():assert json.loads(path.read_text())==result;checked+=1
   else:path.write_text(json.dumps(result,separators=(',',':'))+'\n')
  record=dict(status='passed',sourceURL=URL,sourceETag=ETAG,sourceBytes=SIZE,nodes=len(CACHE),indexBytes=sum(len(x) for x in CACHE.values()),hdf5Version=h5py.version.hdf5_version,chunkCounts={k:len(v['values']) for k,v in maps.items()},allRuns=3456,earlierCoordinateMapsMatched=checked,seconds=time.time()-before,requests=RECEIPTS,otherHDF5Reads=f.requests)
  (out/'chunk-index-verification.json').write_text(json.dumps(record,indent=2)+'\n');print('ALL CHUNK MAPS COMPLETE',3456,flush=True)
 finally:(out/'chunk-index-ranges.json').write_text(json.dumps(RECEIPTS,indent=2)+'\n')
