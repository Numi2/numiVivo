#!/usr/bin/env python3
"""Run frozen native null fits and replay; retain source-deduplicated archives."""
import argparse,gzip,hashlib,json,os,shutil,subprocess,time
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--hdf5',type=Path,required=True)
a=p.parse_args();assert sha(a.binary)=='4eee4da0f3cdfef0ea47d341a6ce75d27042d2ff8e017b57be3f5e9602bb1884'
env={**os.environ,'NUMIVIVO_HDF5_LIBRARY':str(a.hdf5)};runs=[]
def run(args,log):
 start=time.monotonic()
 with log.open('w') as f:r=subprocess.run([str(a.binary)]+list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,env=env)
 record=dict(command=list(map(str,args)),seconds=time.monotonic()-start,exitCode=r.returncode,log=str(log))
 runs.append(record);write(a.root/'native-runs.json',dict(binarySHA256=sha(a.binary),hdf5SHA256=sha(a.hdf5),runs=runs))
 print(json.dumps(record),flush=True);return r.returncode
for study in ['kang','hagai']:
 d=a.root/study;source=d/'original.h5ad';wire=json.loads((d/'native-input.json').read_text())
 assert sha(source)==wire['sourceSHA256']
 annotated=d/'annotated.h5ad'
 assert run(['singlecell-h5ad-annotate',source,'--plan',d/'annotation-plan.json','--output',annotated],d/'annotation.log')==0
 source_hash=sha(annotated)
 for seed in range(1,11):
  s=d/str(seed);s.mkdir(exist_ok=False)
  combined=wire['plans'][str(seed)]
  for mode_index,mode in enumerate(['default','active']):
   assert shutil.disk_usage(d).free>2*annotated.stat().st_size+100_000_000,'Insufficient space for next native snapshot/replay'
   plan=dict(combined,contrasts=[combined['contrasts'][mode_index]],cellSelection=dict(source=dict(bytes=list(bytes.fromhex(source_hash))),observationIndices=wire['sourceRows'],provenance='Frozen untreated sham cohort, protocol SHA256='+wire['protocolSHA256']))
   path=s/(mode+'-plan.json');write(path,plan);out=s/mode
   code=run(['singlecell-h5ad-pseudobulk',annotated,'--plan',path,'--output',out],s/(mode+'-publish.log'))
   if code:assert not out.exists();continue
   verified=run(['singlecell-h5ad-pseudobulk-verify',out],s/(mode+'-verify.log'))==0
   assert verified,'Do not archive a failed replay as verified'
   assert sha(out/'original.h5ad')==source_hash
   entries=[]
   for name in ['plan.json','report.json','receipt.json']:
    f=out/name;raw=f.read_bytes();compressed=gzip.compress(raw,mtime=0);target=out/(name+'.gz');target.write_bytes(compressed)
    assert gzip.decompress(target.read_bytes())==raw
    entries.append(dict(path=target.name,sha256=sha(target),logicalPath=name,logicalSHA256=hashlib.sha256(raw).hexdigest()))
   handles=subprocess.run(['/usr/sbin/lsof','-F','p','--',str(out/'original.h5ad')],capture_output=True,text=True)
   assert handles.returncode==1 and not handles.stdout and not handles.stderr
   write(out/'archive.json',dict(status='native-publication-and-replay-verified-before-archival',source='../../annotated.h5ad',sourceSHA256=source_hash,entries=entries,
    completeSourceRetained=True,noOpenSourceHandles=True,restoreInstruction='Use restore_bundle.py to materialize a normal verifiable bundle in a new directory.'))
   for name in ['original.h5ad','plan.json','report.json','receipt.json']:(out/name).unlink()
write(a.root/'native-complete.json',dict(status='finished-all-declared-attempts',runs=len(runs),failedCommands=sum(x['exitCode']!=0 for x in runs)))
