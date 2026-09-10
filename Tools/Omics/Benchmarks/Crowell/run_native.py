#!/usr/bin/env python3
"""Run every predeclared population; retain unavailable replication and fit failures."""
import argparse,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','binary','hdf5']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();protocol=json.loads((a.root/'protocol.json').read_text());source=a.root/'counts.h5ad'
assert sha(a.binary)=='2d9e524cbb1cdd552fc9f52a1c5237cc904f9faec6eb6bc30b7f49d479ee5433';assert sha(source)==protocol['sourceH5ADSHA256']
env={**os.environ,'NUMIVIVO_HDF5_LIBRARY':str(a.hdf5)};runs=[]
def write(p,d):p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
def run(case,stage,args):
 log=a.root/case['id']/(stage+'.log');start=time.monotonic()
 with log.open('w') as f:r=subprocess.run(['/usr/bin/time','-l',str(a.binary)]+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,env=env)
 rec=dict(case=case['id'],cellGroup=case['cellGroup'],stage=stage,exitCode=r.returncode,seconds=time.monotonic()-start,expectedReplicationGate=case['expectedReplicationGate'],arguments=list(map(str,args)))
 runs.append(rec);write(a.root/'native-runs.json',dict(binarySHA256=sha(a.binary),hdf5SHA256=sha(a.hdf5),sourceSHA256=sha(source),runs=runs));print(json.dumps(rec),flush=True);return r.returncode
for case in protocol['cases']:
 d=a.root/case['id'];out=d/'native';assert sha(d/'plan.json')==case['planSHA256']
 code=run(case,'publish',['singlecell-h5ad-pseudobulk',source,'--plan',d/'plan.json','--output',out])
 if code:
  assert not out.exists();continue
 assert case['expectedReplicationGate']=='ready','Unexpected admission of an insufficient animal cohort'
 assert run(case,'replay',['singlecell-h5ad-pseudobulk-verify',out])==0
 assert sha(out/'original.h5ad')==protocol['sourceH5ADSHA256']
 entries=[]
 for name in ['plan.json','report.json','receipt.json']:
  raw=(out/name).read_bytes();target=out/(name+'.gz');target.write_bytes(gzip.compress(raw,mtime=0));assert gzip.decompress(target.read_bytes())==raw
  entries.append(dict(path=target.name,sha256=sha(target),logicalPath=name,logicalSHA256=hashlib.sha256(raw).hexdigest()))
 handles=subprocess.run(['lsof',str(out/'original.h5ad')],capture_output=True);assert handles.returncode==1 and not handles.stdout
 write(out/'archive.json',dict(status='published-and-replay-verified-before-archival',source=str(source),sourceSHA256=sha(source),entries=entries,completeSourceRetained=True,restoreTool='Tools/Omics/NegativeBinomial/EffectShrinkage/restore_bundle.py --source counts.h5ad'))
 for name in ['original.h5ad','plan.json','report.json','receipt.json']:(out/name).unlink()
unexpected=[r for r in runs if r['exitCode'] and r['expectedReplicationGate']=='ready']
write(a.root/'native-complete.json',dict(status='all-declared-populations-attempted',commands=len(runs),successfulPublications=sum(r['stage']=='publish' and r['exitCode']==0 for r in runs),successfulReplays=sum(r['stage']=='replay' and r['exitCode']==0 for r in runs),replicationGateFailures=[r for r in runs if r['exitCode'] and r['expectedReplicationGate']=='insufficient'],unexpectedFailures=unexpected))
if unexpected:raise SystemExit(1)
