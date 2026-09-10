#!/usr/bin/env python3
import argparse,concurrent.futures,gzip,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--ql-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--jobs',type=int,default=2);a=p.parse_args();assert 1<=a.jobs<=3
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();script=Path(__file__).with_suffix('.R');protocol=json.loads((a.ql_root/'protocol.json').read_text())
def execute(c):
 source=a.ql_root/c['id'];out=a.out/c['id'];out.mkdir(parents=True,exist_ok=True)
 inputSHA=sha(source/'input.json.gz');priorSHA=sha(source/'reference.json.gz');dest=out/'reference.json.gz'
 if (out/'reference-receipt.json').exists():
  r=json.loads((out/'reference-receipt.json').read_text());assert r['inputSHA256']==inputSHA and r['priorSHA256']==priorSHA and r['scriptSHA256']==sha(script) and r['outputSHA256']==sha(dest);return r
 with (out/'reference.log').open('wb') as log:
  run=subprocess.run(['/opt/homebrew/bin/Rscript',str(script),str(source/'input.json.gz'),str(source/'reference.json.gz'),str(out/'reference.json')],stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'R_LIBS_USER':'/Users/home/numivivo-r-library-20260909','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
 r=dict(case=c['id'],exitCode=run.returncode,inputSHA256=inputSHA,priorSHA256=priorSHA,scriptSHA256=sha(script),protocolSHA256=sha(script.with_name('PROTOCOL.md')))
 if (out/'reference.json').exists():
  raw=(out/'reference.json').read_bytes();dest.write_bytes(gzip.compress(raw,mtime=0));assert gzip.decompress(dest.read_bytes())==raw;(out/'reference.json').unlink();r.update(outputSHA256=sha(dest),logicalSHA256=hashlib.sha256(raw).hexdigest())
 (out/'reference-receipt.json').write_text(json.dumps(r,sort_keys=True,indent=2)+'\n');print(json.dumps(r),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as executor:results=list(executor.map(execute,protocol['cases']))
(a.out/'reference-complete.json').write_text(json.dumps(results,sort_keys=True,indent=2)+'\n')
