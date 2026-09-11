#!/usr/bin/env python3
"""Run frozen native intervals and preserve verified compressed bundles to bound disk use."""
import argparse,hashlib,json,os,shutil,subprocess,tarfile,tempfile,time
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);(a.out/'plans').mkdir();(a.out/'bundles').mkdir()
 prior=json.loads((a.inputs/'input-freeze.json').read_text())
 for name,digest in prior['files'].items():assert sha(a.inputs/name)==digest,name
 folds=json.loads((a.inputs/'folds.json').read_text());assert len(folds)==26
 for f in folds:
  plan=json.loads((a.inputs/f['id']/'training.json').read_text());assert 'donorResponseIntervalCoverage' not in plan;plan['donorResponseIntervalCoverage']=0.95;write(a.out/'plans'/(f['id']+'.json'),plan)
 write(a.out/'input-freeze.json',dict(schemaVersion=1,priorInputFreezeSHA256=sha(a.inputs/'input-freeze.json'),files={str(p.relative_to(a.out)):sha(p) for p in sorted((a.out/'plans').glob('*.json'))},protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),runnerSHA256=sha(__file__),binarySHA256=sha(a.binary),fitStarted=False))
 write(a.out/'execution.json',dict(startedUnix=time.time(),binary=str(a.binary),binarySHA256=sha(a.binary),sourceInputs=str(a.inputs),inputFreezeSHA256=sha(a.out/'input-freeze.json')))
 commands=[];archives=[]
 def run(args,log):
  start=time.time();r=subprocess.run([str(a.binary.resolve()),*map(str,args)],capture_output=True,text=True);log.write_text(r.stdout+r.stderr);commands.append(dict(arguments=list(map(str,args)),exitCode=r.returncode,elapsedSeconds=time.time()-start));write(a.out/'commands.json',commands);assert r.returncode==0,(args,r.stdout,r.stderr)
 for index,f in enumerate(folds):
  root=Path(tempfile.mkdtemp(prefix='scratch-',dir=a.out));src=a.inputs/f['id']
  # No automatic deletion on error: retain failed native evidence and diagnose it.
  run(['singlecell-perturbation-fit',src/'training.h5ad','--plan',a.out/'plans'/(f['id']+'.json'),'--output',root/'model'],root/'fit.log')
  run(['singlecell-perturbation-verify',root/'model'],root/'verify-model.log')
  args=['singlecell-perturbation-predict',src/'query.h5ad','--plan',src/'query.json','--reference',root/'model','--output']
  run(args+[root/'prediction'],root/'predict.log');run(['singlecell-perturbation-prediction-verify',root/'prediction'],root/'verify-prediction.log')
  if index==0:
   run(args+[root/'repeat'],root/'repeat.log');assert (root/'prediction/report.json').read_bytes()==(root/'repeat/report.json').read_bytes()
  files=sorted(p for p in root.rglob('*') if p.is_file());assert all(not p.is_symlink() for p in root.rglob('*'))
  members={str(p.relative_to(root)):dict(bytes=p.stat().st_size,SHA256=sha(p)) for p in files};archive=a.out/'bundles'/(f['id']+'.tar.gz')
  with tarfile.open(archive,'w:gz',compresslevel=6) as t:
   for p in files:t.add(p,arcname=str(p.relative_to(root)),recursive=False)
  with tarfile.open(archive,'r:gz') as t:
   items=t.getmembers();assert len(items)==len(members) and {x.name for x in items}==set(members)
   for item in items:
    assert item.isfile();b=t.extractfile(item).read();assert len(b)==members[item.name]['bytes'] and hashlib.sha256(b).hexdigest()==members[item.name]['SHA256']
  opened=subprocess.run(['/usr/sbin/lsof','+D',str(root)],capture_output=True,text=True);assert opened.returncode==1 and not opened.stdout.strip() and not opened.stderr.strip()
  archives.append(dict(fold=f['id'],path=str(archive.relative_to(a.out)),bytes=archive.stat().st_size,SHA256=sha(archive),members=members,uncompressedBytes=sum(x['bytes'] for x in members.values()),scratchPath=str(root),archiveVerified=True))
  write(a.out/'archive-progress.json',archives)
  # Only this finished invocation's scratch is removed, after complete archive
  # byte verification and a targeted open-handle check. Source inputs are retained.
  shutil.rmtree(root);print(json.dumps(dict(completed=index+1,fold=f['id'],archiveBytes=archive.stat().st_size)),flush=True)
 write(a.out/'prediction-freeze.json',dict(schemaVersion=1,inputFreezeSHA256=sha(a.out/'input-freeze.json'),folds=archives,commands=len(commands),scoringStarted=False,completedUnix=time.time()))
 print(json.dumps(dict(status='completed',folds=len(archives),commands=len(commands),archiveBytes=sum(x['bytes'] for x in archives),originalBytes=sum(x['uncompressedBytes'] for x in archives),freezeSHA256=sha(a.out/'prediction-freeze.json'))),flush=True)
if __name__=='__main__':main()
