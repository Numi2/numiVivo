from pathlib import Path
import json,hashlib,subprocess,os,time,re
r=Path(__file__).parent;parent=Path('/Users/n/numivivo-metal-kang-20260912');source=parent/'pca/scores.bin';metadata=parent/'pca/metadata.json';python=r/'env/bin/python'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
original=json.loads((parent/'protocol.json').read_text());assert sha(source)==original['originalPCASHA256'];assert sha(metadata)==original['originalMetadataSHA256']
protocol=dict(order=['cpu','scanpy','metal','scanpy','metal','cpu','metal','cpu','scanpy'],cells=24673,components=20,neighborsIncludingSelf=20,cpuWorkers=1,sourceSHA256=sha(source),metadataSHA256=sha(metadata),nativeBinarySHA256=sha(r/'graph'),driverSHA256=sha(Path(__file__)),scanpyDriverSHA256=sha(r/'scanpy_graph.py'),nativeDriverSHA256=sha(r/'Graph.swift'),scope='Same frozen PCA graph construction only: native CPU FP64, native Metal FP32 and exact sklearn-backed Scanpy UMAP fuzzy graph. Source reads included; output serialization, imports and startup excluded from graphStageSeconds but retained in process wall time. No whole-CLI speedup inference.',createdUnix=time.time())
(r/'protocol.json').write_text(json.dumps(protocol,indent=2));results=[]
for i,backend in enumerate(protocol['order']):
 name=str(i)+'-'+backend;dest=r/name;args=[str(python),str(r/'scanpy_graph.py'),str(source),str(dest)] if backend=='scanpy' else [str(r/'graph'),backend,str(source),str(metadata),str(dest)]
 start=time.time()
 with (r/(name+'.log')).open('w') as f:code=subprocess.call(['/usr/bin/time','-l']+args,stdout=f,stderr=subprocess.STDOUT)
 elapsed=time.time()-start;assert code==0,name
 timing=json.loads((dest/'timing.json').read_text());log=(r/(name+'.log')).read_text();results.append(dict(run=i,backend=backend,exitCode=code,graphStageSeconds=timing['graphStageSeconds'],processWallSeconds=elapsed,peakRSSBytes=int(re.search(r'(\d+)\s+maximum resident set size',log)[1])));(r/'status.json').write_text(json.dumps(results,indent=2))
(r/'terminal.json').write_text(json.dumps(dict(results=results,completedUnix=time.time()),indent=2))
