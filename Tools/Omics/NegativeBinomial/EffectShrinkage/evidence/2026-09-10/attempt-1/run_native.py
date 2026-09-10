#!/usr/bin/env python3
"""Publish and replay the four frozen treatment contrasts, then losslessly archive."""
import argparse,gzip,hashlib,json,os,subprocess,time
from pathlib import Path

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d): p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','binary','hdf5','sources']: p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();runs=[];binary_hash=sha(a.binary)
protocol=json.loads((a.root/'protocol.json').read_text())
env={**os.environ,'NUMIVIVO_HDF5_LIBRARY':str(a.hdf5)}
for rec in protocol['records']:
 d=a.root/(rec['study']+'-'+rec['mode']);source=a.sources/rec['study']/'original.h5ad';out=d/'native'
 assert sha(source)==rec['sourceSHA256']
 for name,args in [('publish',['singlecell-h5ad-pseudobulk',source,'--plan',d/'plan.json','--output',out]),('replay',['singlecell-h5ad-pseudobulk-verify',out])]:
  assert sha(a.binary)==binary_hash
  start=time.monotonic()
  with (d/(name+'.log')).open('w') as f:
   r=subprocess.run(['/usr/bin/time','-l',str(a.binary)]+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,env=env)
  row=dict(study=rec['study'],mode=rec['mode'],stage=name,seconds=time.monotonic()-start,exitCode=r.returncode,command=list(map(str,args)))
  runs.append(row);write(a.root/'native-runs.json',dict(binarySHA256=binary_hash,hdf5SHA256=sha(a.hdf5),protocolSHA256=sha(a.root/'protocol.json'),runs=runs));print(json.dumps(row),flush=True)
  assert r.returncode==0,'Retain failed attempt; do not archive as successful'
 assert sha(out/'original.h5ad')==rec['sourceSHA256']
 entries=[]
 for name in ['plan.json','report.json','receipt.json']:
  raw=(out/name).read_bytes();target=out/(name+'.gz');target.write_bytes(gzip.compress(raw,mtime=0))
  assert gzip.decompress(target.read_bytes())==raw
  entries.append(dict(path=target.name,sha256=sha(target),logicalPath=name,logicalSHA256=hashlib.sha256(raw).hexdigest()))
 check=subprocess.run(['/usr/sbin/lsof',str(out/'original.h5ad')],capture_output=True);assert check.returncode==1 and not check.stdout
 write(out/'archive.json',dict(status='published-and-replay-verified-before-archival',source=str(source),sourceSHA256=sha(source),entries=entries,restoreTool='../NullBenchmark/restore_bundle.py',completeSourceRetained=True))
 for name in ['original.h5ad','plan.json','report.json','receipt.json']:(out/name).unlink()
write(a.root/'native-complete.json',dict(status='four-treatment-contrasts-published-and-replayed',commands=len(runs),binarySHA256=binary_hash))
