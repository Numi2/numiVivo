#!/usr/bin/env python3
"""Run all frozen cases, retaining every reference stage and terminal result."""
import argparse,concurrent.futures,gzip,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--r-library',type=Path,required=True);p.add_argument('--jobs',type=int,default=2);a=p.parse_args()
script=Path(__file__).with_name('reference.R');protocol=json.loads((a.root/'protocol.json').read_text());assert 1<=a.jobs<=4
assert hashlib.sha256(script.with_name('PROTOCOL.md').read_bytes()).hexdigest()==protocol['protocolSHA256']
def run(c):
 d=a.root/c['id'];out=d/'reference.json';assert not out.exists() and not (d/'reference.json.gz').exists()
 with (d/'reference.log').open('w') as log:r=subprocess.run(['/opt/homebrew/bin/Rscript',str(script),str(d/'input.json.gz'),str(out)],stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'R_LIBS_USER':str(a.r_library),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
 rec=dict(case=c['id'],exitCode=r.returncode,scriptSHA256=hashlib.sha256(script.read_bytes()).hexdigest())
 if out.exists():
  raw=out.read_bytes();packed=gzip.compress(raw,mtime=0);(d/'reference.json.gz').write_bytes(packed);assert gzip.decompress(packed)==raw;out.unlink();rec.update(outputSHA256=hashlib.sha256(packed).hexdigest(),logicalSHA256=hashlib.sha256(raw).hexdigest())
 (d/'run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:results=list(ex.map(run,protocol['cases']))
(a.root/'complete.json').write_text(json.dumps(results,indent=2)+'\n')
