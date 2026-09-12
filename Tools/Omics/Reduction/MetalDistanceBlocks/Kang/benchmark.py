from pathlib import Path
import json,hashlib,subprocess,os,time,re,statistics
r=Path(__file__).parent;binary=Path('/Users/n/numivivo-native-metal-knn-20260912/build/numivivo-omics')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
p=json.loads((r/'protocol.json').read_text());assert sha(binary)==p['binarySHA256'];out=r/'benchmark';out.mkdir()
protocol=dict(order=['cpu','metal','metal','cpu','cpu','metal'],binarySHA256=sha(binary),planSHA256={m:sha(r/(m+'-plan.json')) for m in ['cpu','metal']},inputReceiptSHA256=sha(r/'pca/receipt.json'),scope='Complete native graph-publication CLI, including parent PCA reconstruction; warm exercised caches; CPU FP64 versus Metal FP32; no Scanpy performance claim',createdUnix=time.time());(out/'protocol.json').write_text(json.dumps(protocol,indent=2));results=[]
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-integration-reference-py/lib/python3.13/site-packages/h5py/.dylibs/libhdf5.320.0.0.dylib')
for i,mode in enumerate(protocol['order']):
 name=str(i)+'-'+mode;dest=out/name;start=time.time()
 with (out/(name+'.log')).open('w') as f:code=subprocess.call(['/usr/bin/time','-l',str(binary),'singlecell-pca-neighbors',str(r/'pca'),'--plan',str(r/(mode+'-plan.json')),'--output',str(dest)],stdout=f,stderr=subprocess.STDOUT,env=env)
 elapsed=time.time()-start;assert code==0
 files={n:sha(dest/n) for n in ['neighbors.bin','edges.bin','offsets.bin','bandwidths.bin','graph.json','plan.json','execution.json','receipt.json']}
 assert all(h==sha(r/mode/n) for n,h in files.items())
 log=(out/(name+'.log')).read_text();rss=int(re.search(r'(\d+)\s+maximum resident set size',log)[1]);results.append(dict(run=i,backend=mode,exitCode=code,processWallSeconds=elapsed,peakRSSBytes=rss,allEightFilesMatch=True,hashes=files));(out/'status.json').write_text(json.dumps(results,indent=2))
summary=dict(results=results,medianSeconds={m:statistics.median(x['processWallSeconds'] for x in results if x['backend']==m) for m in ['cpu','metal']},meanSeconds={m:statistics.mean(x['processWallSeconds'] for x in results if x['backend']==m) for m in ['cpu','metal']},completedUnix=time.time());(out/'terminal.json').write_text(json.dumps(summary,indent=2));print(summary['medianSeconds'])
