#!/usr/bin/env python3
"""Stream native measurements from immutable remote source archives."""
import argparse,concurrent.futures,gzip,hashlib,json,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
p.add_argument('--binary',required=True);p.add_argument('--jobs',type=int,default=4);a=p.parse_args()
protocol=json.loads((a.root/'protocol.json').read_text());host=protocol['host']
assert hashlib.sha256(Path(__file__).with_name('PROTOCOL.md').read_bytes()).hexdigest()==protocol['protocolSHA256']
def remote_sha(path):return subprocess.check_output(['ssh',host,'shasum','-a','256',path],text=True).split()[0]
binary_sha=remote_sha(a.binary)
def run(rec):
 d=a.root/rec['id'];d.mkdir(exist_ok=False);assert remote_sha(rec['source'])==rec['sourceSHA256']
 assert remote_sha(a.binary)==binary_sha
 command=['ssh',host,'/usr/bin/time','-l',a.binary,rec['source']];start=time.monotonic()
 with (d/'native.log').open('wb') as log:r=subprocess.run(command,stdout=subprocess.PIPE,stderr=log)
 result=dict(case=rec['id'],command=command,exitCode=r.returncode,seconds=time.monotonic()-start,binarySHA256=binary_sha,sourceSHA256=rec['sourceSHA256'])
 if r.returncode==0:
  json.loads(r.stdout);raw=gzip.compress(r.stdout,mtime=0);(d/'native.json.gz').write_bytes(raw)
  result.update(outputSHA256=hashlib.sha256(raw).hexdigest(),logicalSHA256=hashlib.sha256(r.stdout).hexdigest())
 (d/'run.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(dict(case=rec['id'],exitCode=r.returncode)),flush=True)
 return result
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:results=list(pool.map(run,protocol['cases']))
(a.root/'native-complete.json').write_text(json.dumps(dict(binarySHA256=binary_sha,runs=results),sort_keys=True,indent=2)+'\n')
assert all(r['exitCode']==0 for r in results)
