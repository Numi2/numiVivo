from pathlib import Path
import json,subprocess,time,hashlib
s=Path(__file__).parent;freeze=json.loads((s/'input-freeze.json').read_text());out=s/'native';out.mkdir(exist_ok=False);record={'scoringStarted':False,'inputFreezeSHA256':hashlib.sha256((s/'input-freeze.json').read_bytes()).hexdigest(),'binarySHA256':hashlib.sha256((s/'study-shrinkage').read_bytes()).hexdigest(),'folds':{}}
for name,digest in freeze['inputs'].items():
 p=s/'inputs'/name;assert hashlib.sha256(p.read_bytes()).hexdigest()==digest;dest=out/name;start=time.time();subprocess.run([str(s/'study-shrinkage'),str(p),str(dest)],check=True);record['folds'][p.stem]={'inputSHA256':digest,'outputSHA256':hashlib.sha256(dest.read_bytes()).hexdigest(),'seconds':time.time()-start};print(p.stem,flush=True)
record['status']='complete-three-study-native-predictions';(s/'prediction-freeze.json').write_text(json.dumps(record,indent=2)+'\n')
