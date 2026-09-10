#!/usr/bin/env python3
"""Publish, replay and losslessly archive each frozen empirical-prior contrast."""
import argparse,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','binary','hdf5']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();runs=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def write(p,d):p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
binary_hash=sha(a.binary);protocol=json.loads((a.root/'protocol.json').read_text());env={**os.environ,'NUMIVIVO_HDF5_LIBRARY':str(a.hdf5)}
for rec in protocol['records']:
 d=a.root/rec['id'];out=d/'native';source=Path(rec['source']);assert sha(source)==rec['sourceSHA256'];assert sha(d/'plan.json')==rec['planSHA256']
 for stage,args in [('publish',['singlecell-h5ad-pseudobulk',source,'--plan',d/'plan.json','--output',out]),('replay',['singlecell-h5ad-pseudobulk-verify',out])]:
  assert sha(a.binary)==binary_hash;start=time.monotonic()
  with (d/(stage+'.log')).open('w') as f:r=subprocess.run(['/usr/bin/time','-l',str(a.binary)]+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,env=env)
  row=dict(case=rec['id'],stage=stage,exitCode=r.returncode,seconds=time.monotonic()-start,command=list(map(str,args)));runs.append(row)
  write(a.root/'native-runs.json',dict(binarySHA256=binary_hash,hdf5SHA256=sha(a.hdf5),runs=runs));print(json.dumps(row),flush=True);assert r.returncode==0
 assert sha(out/'original.h5ad')==rec['sourceSHA256'];entries=[]
 for name in ['plan.json','report.json','receipt.json']:
  raw=(out/name).read_bytes();target=out/(name+'.gz');target.write_bytes(gzip.compress(raw,mtime=0));assert gzip.decompress(target.read_bytes())==raw
  entries.append(dict(path=target.name,sha256=sha(target),logicalPath=name,logicalSHA256=hashlib.sha256(raw).hexdigest()))
 handles=subprocess.run(['lsof','+D',str(out)],capture_output=True);assert handles.returncode==1 and not handles.stdout
 write(out/'archive.json',dict(status='published-and-replay-verified-before-archival',source=str(source),sourceSHA256=sha(source),entries=entries,completeSourceRetained=True,restoreTool='Tools/Omics/NegativeBinomial/EffectShrinkage/restore_bundle.py'))
 for name in ['original.h5ad','plan.json','report.json','receipt.json']:(out/name).unlink()
write(a.root/'native-complete.json',dict(cases=len(protocol['records']),commands=len(runs),status='all-frozen-contrasts-published-and-replayed',binarySHA256=binary_hash))
