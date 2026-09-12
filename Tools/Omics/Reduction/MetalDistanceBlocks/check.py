from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import numpy as np
from scipy.spatial.distance import cdist
import json,hashlib,re,statistics
r=Path(__file__).parent;p=json.loads((r/'protocol.json').read_text());terminal=json.loads((r/'terminal.json').read_text());assert len(terminal['results'])==6
assert all(x['exitCode']==0 for x in terminal['results'])
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
source=Path(p['source']);assert sha(source)==p['sourceSHA256']
n,d,k=p['cells'],p['dimensions'],p['neighborsIncludingSelf']
a=np.fromfile(source,dtype=[('row','<u4'),('col','<u4'),('value','<f8')]);assert np.array_equal(a['row'],np.repeat(np.arange(n),d));assert np.array_equal(a['col'],np.tile(np.arange(d),n));x=a['value'].reshape(n,d)
outputs={};timings=[]
for entry in terminal['results']:
 i=entry['run'];mode=entry['backend'];f=r/(str(i)+'-'+mode);ids=np.fromfile(f/'indices.bin','<u4').reshape(n,k);v=np.fromfile(f/'squared-distances.bin','<f4').reshape(n,k)
 assert np.array_equal(ids[:,0],np.arange(n)) and np.all(v[:,0]==0) and np.all(ids<n) and np.all(np.isfinite(v))
 assert all(len(set(row))==k for row in ids)
 if mode in outputs:assert np.array_equal(ids,outputs[mode][0]) and np.array_equal(v,outputs[mode][1])
 else:outputs[mode]=(ids,v)
 report=json.loads((f/'report.json').read_text());log=(r/(str(i)+'-'+mode+'.log')).read_text();report['peakRSSBytes']=int(re.search(r'(\d+)\s+maximum resident set size',log)[1]);report['processWallSeconds']=entry['wallSeconds'];timings.append(report)
comparison={mode:dict(changedRows=0,missingNeighbors=0,maximumAbsoluteSquaredDistanceError=0.,maximumScaledSquaredDistanceError=0.,differences=[]) for mode in outputs}
for start in range(0,n,128):
 stop=min(n,start+128);dist=cdist(x[start:stop],x,'sqeuclidean');dist[np.arange(stop-start),np.arange(start,stop)]=np.inf
 thresholds=np.partition(dist,k-2,axis=1)[:,k-2]
 for local,row in enumerate(range(start,stop)):
  candidates=np.flatnonzero(dist[local]<=thresholds[local]);order=np.lexsort((candidates,dist[local,candidates]));expected=candidates[order[:k-1]]
  for mode,(ids,v) in outputs.items():
   z=comparison[mode];got=ids[row,1:];missing=sorted(set(expected)-set(got));extra=sorted(set(got)-set(expected))
   if not np.array_equal(got,expected):z['changedRows']+=1;z['differences'].append(dict(row=row,missing=list(map(int,missing)),extra=list(map(int,extra)),orderOnly=not missing))
   z['missingNeighbors']+=len(missing)
   truth=dist[local,got];err=abs(v[row,1:].astype(float)-truth)
   z['maximumAbsoluteSquaredDistanceError']=max(z['maximumAbsoluteSquaredDistanceError'],float(err.max()))
   z['maximumScaledSquaredDistanceError']=max(z['maximumScaledSquaredDistanceError'],float((err/np.maximum(1,truth)).max()))
for z in comparison.values():z['neighborRecall']=1-z['missingNeighbors']/(n*(k-1));assert z['maximumScaledSquaredDistanceError']<1e-5
result=dict(scope=p['scope'],allCells=n,allDirectedDistances=n*(n-1),reference='SciPy float64 sqeuclidean; complete candidates and row-ID tie breaking',comparison=comparison,cpuMetalIndicesExact=bool(np.array_equal(outputs['cpu'][0],outputs['metal'][0])),cpuMetalDistancesExact=bool(np.array_equal(outputs['cpu'][1],outputs['metal'][1])),timings=timings,medianProcessSeconds={m:statistics.median(t['processWallSeconds'] for t in timings if t['backend']==m) for m in outputs},sourceSHA256=p['sourceSHA256'])
(r/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps({key:value for key,value in result.items() if key not in ['timings','comparison']}));print(json.dumps({m:{k:v for k,v in z.items() if k!='differences'} for m,z in comparison.items()}))
