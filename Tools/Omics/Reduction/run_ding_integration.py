#!/usr/bin/env python3
"""Run full Ding cohorts remotely, verify and retain each artifact before disposal."""
import argparse,hashlib,json,re,shlex,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--host',required=True)
for name in ['binary','hdf5','remote-root','prepared','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[];artifacts=[]
protocol=json.loads((a.prepared/'protocol.json').read_text());source=json.loads((a.prepared/'checks.json').read_text())
def save(name,v):(a.out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def remote(args):return subprocess.run(['ssh',a.host,shlex.join(map(str,args))],capture_output=True,text=True)
def run(label,args,ok=True):
 command=['env','NUMIVIVO_HDF5_LIBRARY='+str(a.hdf5),'/usr/bin/time','-l',a.binary,*args]
 r=remote(command);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
 rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);seconds=re.search(r'([\d.]+) real',r.stderr)
 commands.append(dict(label=label,command=list(map(str,command)),exitCode=r.returncode,expectedSuccess=ok,seconds=float(seconds[1]) if seconds else None,maximumResidentBytes=int(rss[1]) if rss else None));save('commands.json',commands)
 assert r.returncode==(0 if ok else 65),(label,r.stderr);print(label+' passed',flush=True);return r
def transfer(remote_path,local_path):
 local_path.mkdir(parents=True,exist_ok=False)
 subprocess.run(['rsync','-az','--exclude','*.h5ad',a.host+':'+str(remote_path)+'/',str(local_path)+'/'],check=True)
 script='''import hashlib,json,sys
from pathlib import Path
r=Path(sys.argv[1])
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
print(json.dumps([dict(path=p.relative_to(r).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(r.rglob('*')) if p.is_file()]))
'''
 result=remote(['python3','-c',script,remote_path]);assert result.returncode==0,result.stderr
 entries=json.loads(result.stdout)
 for e in entries:
  if e['path'].endswith('.h5ad'):assert e['sha256']==source['preparedSHA256'];continue
  path=local_path/e['path'];assert path.stat().st_size==e['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest()==e['sha256'],path
 artifacts.append(dict(remotePath=str(remote_path),localPath=str(local_path),files=entries,allTransferredBytesExact=True));save('artifacts.json',artifacts)
 return entries
binary=remote(['shasum','-a','256',a.binary]);assert binary.returncode==0
binary_hash=binary.stdout.split()[0];assert binary_hash=='a56a056349aac2d2463ad56b5a8260257186cf6388a021f77d5bb0bc38566c8d'
result=remote(['mkdir',a.remote_root]);assert result.returncode==0,result.stderr
subprocess.run(['rsync','-az',str(a.prepared)+'/',a.host+':'+str(a.remote_root/'prepared')+'/'],check=True)
remote_source=remote(['shasum','-a','256',a.remote_root/'prepared/ding.h5ad']);assert remote_source.returncode==0 and remote_source.stdout.split()[0]==source['preparedSHA256']
run('fit',['singlecell-h5ad-pca',a.remote_root/'prepared/ding.h5ad','--plan',a.remote_root/'prepared/fit.json','--output',a.remote_root/'pca'])
transfer(a.remote_root/'pca',a.out/'pca')
for mode,extra in protocol['modes'].items():
 for seed in protocol['seeds']:
  label=mode+'-'+str(seed);options=dict(protocol['integration'],**extra,seed=seed)
  plan=dict(schemaVersion=1,inputKind='fitted',integration=options);local_plan=a.out/(label+'-plan.json');local_plan.write_text(json.dumps(plan,indent=2)+'\n')
  remote_plan=a.remote_root/(label+'-plan.json');subprocess.run(['rsync',str(local_plan),a.host+':'+str(remote_plan)],check=True)
  target=a.remote_root/label
  run(label,['singlecell-pca-integrate',a.remote_root/'pca','--plan',remote_plan,'--output',target]);run(label+'-verify',['singlecell-pca-integrate-verify',target])
  entries=transfer(target,a.out/label)
  # Only this newly created, fully verified run is disposable. Check its exact
  # inventory and absence of open handles; preserve the canonical source/PCA.
  expected=a.out/(label+'-inventory.json');expected.write_text(json.dumps(entries)+'\n');subprocess.run(['rsync',str(expected),a.host+':'+str(a.remote_root/'disposable-inventory.json')],check=True)
  script='''import hashlib,json,shutil,subprocess,sys
from pathlib import Path
r=Path(sys.argv[1]);root=Path(sys.argv[2]);assert r.parent==root and r.name in ['fixed-7','fixed-19','fixed-41','adaptive-7','adaptive-19','adaptive-41'] and not r.is_symlink()
expected={e['path']:e for e in json.loads((root/'disposable-inventory.json').read_text())};files=[p for p in r.rglob('*') if p.is_file()]
assert {p.relative_to(r).as_posix() for p in files}==set(expected)
handles=subprocess.run(['/usr/sbin/lsof','+D',str(r)],capture_output=True,text=True);assert handles.returncode==1 and not handles.stdout and not handles.stderr
for p in files:
 e=expected[p.relative_to(r).as_posix()];assert not p.is_symlink() and p.stat().st_size==e['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==e['sha256']
shutil.rmtree(r)
print(json.dumps(dict(removedVerifiedDisposableRun=str(r),freeBytes=shutil.disk_usage(root).free)))
'''
  result=remote(['python3','-c',script,target,a.remote_root]);assert result.returncode==0,result.stderr
  (a.out/(label+'-disposal.json')).write_text(result.stdout)
save('checks.json',dict(status='passed',binarySHA256=binary_hash,sourceSHA256=source['preparedSHA256'],commands=len(commands),runs=6,allRunsReplayedAndTransferredExactly=True,qualification='All source UMI cells retained. Native numerical execution is separate from source-label biological acceptance. Each full run was verified and backed up before remote disposal; original source/PCA remain available.'))
