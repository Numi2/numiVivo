#!/usr/bin/env python3
"""Publish, replay and losslessly archive each frozen empirical-prior contrast."""
import argparse,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','binary','hdf5']:p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--resume',action='store_true',help='Resume a confirmed terminal attempt; preserve prior command records and logs')
a=p.parse_args();runs=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def write(p,d):p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
binary_hash=sha(a.binary);protocol=json.loads((a.root/'protocol.json').read_text());env={**os.environ,'NUMIVIVO_HDF5_LIBRARY':str(a.hdf5)}
if a.resume:
 previous=json.loads((a.root/'native-runs.json').read_text());assert previous['binarySHA256']==binary_hash and previous['hdf5SHA256']==sha(a.hdf5);runs=previous['runs']
for rec in protocol['records']:
 d=a.root/rec['id'];out=d/'native';source=Path(rec['source']);assert sha(source)==rec['sourceSHA256'];assert sha(d/'plan.json')==rec['planSHA256']
 if a.resume and (out/'archive.json').exists():
  archive=json.loads((out/'archive.json').read_text());assert archive['status']=='published-and-replay-verified-before-archival'
  assert archive['sourceSHA256']==rec['sourceSHA256']
  assert all(any(r['case']==rec['id'] and r['stage']==s and r['exitCode']==0 for r in runs) for s in ['publish','replay'])
  for e in archive['entries']:
   assert sha(out/e['path'])==e['sha256'];assert hashlib.sha256(gzip.decompress((out/e['path']).read_bytes())).hexdigest()==e['logicalSHA256']
  continue
 for stage,args in [('publish',['singlecell-h5ad-pseudobulk',source,'--plan',d/'plan.json','--output',out]),('replay',['singlecell-h5ad-pseudobulk-verify',out])]:
  if a.resume and stage=='publish' and any(r['case']==rec['id'] and r['stage']=='publish' and r['exitCode']==0 for r in runs):
   assert sha(out/'original.h5ad')==rec['sourceSHA256'];continue
  assert sha(a.binary)==binary_hash;start=time.monotonic()
  attempt=1+sum(r['case']==rec['id'] and r['stage']==stage for r in runs)
  log=d/(stage+'.log' if attempt==1 else stage+'-attempt-'+str(attempt)+'.log');assert not log.exists()
  with log.open('w') as f:r=subprocess.run(['/usr/bin/time','-l',str(a.binary)]+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,env=env)
  row=dict(case=rec['id'],stage=stage,attempt=attempt,logPath=str(log),exitCode=r.returncode,seconds=time.monotonic()-start,command=list(map(str,args)));runs.append(row)
  write(a.root/'native-runs.json',dict(binarySHA256=binary_hash,hdf5SHA256=sha(a.hdf5),runs=runs));print(json.dumps(row),flush=True);assert r.returncode==0
 assert sha(out/'original.h5ad')==rec['sourceSHA256'];entries=[]
 for name in ['plan.json','report.json','receipt.json']:
  raw=(out/name).read_bytes();target=out/(name+'.gz');target.write_bytes(gzip.compress(raw,mtime=0));assert gzip.decompress(target.read_bytes())==raw
  entries.append(dict(path=target.name,sha256=sha(target),logicalPath=name,logicalSHA256=hashlib.sha256(raw).hexdigest()))
 handles=subprocess.run(['lsof','+D',str(out)],capture_output=True);assert handles.returncode==1 and not handles.stdout
 write(out/'archive.json',dict(status='published-and-replay-verified-before-archival',source=str(source),sourceSHA256=sha(source),entries=entries,completeSourceRetained=True,restoreTool='Tools/Omics/NegativeBinomial/EffectShrinkage/restore_bundle.py'))
 for name in ['original.h5ad','plan.json','report.json','receipt.json']:(out/name).unlink()
write(a.root/'native-complete.json',dict(cases=len(protocol['records']),commands=len(runs),status='all-frozen-contrasts-published-and-replayed',binarySHA256=binary_hash,retainedFailedCommands=[r for r in runs if r['exitCode']]))
