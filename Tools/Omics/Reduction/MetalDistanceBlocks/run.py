from pathlib import Path
import json,hashlib,subprocess,time
r=Path(__file__).parent;source=Path('/Users/n/numivivo-native-mnn-20260910/hagai/pca/scores.bin')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
prior=json.loads(Path('/Users/n/numivivo-mnn-cli-lifecycle-20260912/qualification.json').read_text())
assert sha(source)==prior['byteEquivalence']['pca/scores.bin']['prior']
protocol=dict(source=str(source),sourceSHA256=sha(source),binarySHA256=sha(r/'knn'),swiftSHA256=sha(r/'Main.swift'),driverSHA256=sha(Path(__file__)),cells=13863,dimensions=20,neighborsIncludingSelf=20,order=['cpu','metal','metal','cpu','cpu','metal'],scope='Full real-Hagai FP32 distance-block research comparison. Independent complete FP64 neighborhood comparison required. No product or biological promotion.',createdUnix=time.time())
(r/'protocol.json').write_text(json.dumps(protocol,indent=2));results=[]
for i,mode in enumerate(protocol['order']):
 assert sha(source)==protocol['sourceSHA256']
 start=time.time()
 with (r/(str(i)+'-'+mode+'.log')).open('w') as f:
  code=subprocess.call(['/usr/bin/time','-l',str(r/'knn'),mode,str(source),'13863','20',str(r/(str(i)+'-'+mode))],stdout=f,stderr=subprocess.STDOUT)
 results.append(dict(run=i,backend=mode,exitCode=code,wallSeconds=time.time()-start));(r/'status.json').write_text(json.dumps(results,indent=2));assert code==0
(r/'terminal.json').write_text(json.dumps(dict(results=results,completedUnix=time.time()),indent=2))
