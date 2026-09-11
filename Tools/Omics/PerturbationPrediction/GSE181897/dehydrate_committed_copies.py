"""Remove only exact committed evidence working copies; preserve all Git blobs."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, time
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--record',type=Path,required=True);a=p.parse_args()
names=[
 'Tools/Omics/PerturbationPrediction/Norman/evidence/2026-09-09-unseen-targets/predictions.npz',
 'Tools/Omics/PerturbationPrediction/Norman/evidence/2026-09-09-control-descriptors/predictions/predictions.npz',
 'Tools/Omics/PerturbationPrediction/Norman/evidence/2026-09-09-native-composition/prediction/report.json.gz',
 'Tools/Omics/PerturbationPrediction/Norman/evidence/2026-09-10-go-transfer/objects/8a1dae1b9ed0408dfe729e6b2b821e2e62c7c3b6282832d377db53bcd7c114ed.gz',
]
def git(*args):return subprocess.check_output(['git','-C',str(a.repo),*args],text=True).strip()
def digest(f):
 h=hashlib.sha256()
 while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def free():s=os.statvfs(a.repo);return s.f_bavail*s.f_frsize
assert not a.record.exists()
head=git('rev-parse','HEAD');manifest=[]
for name in names:
 target=a.repo/name
 assert target.is_file() and not target.is_symlink()
 entry=git('ls-tree',head,'--',name).split()
 assert entry[0]=='100644' and entry[1]=='blob' and entry[3]==name
 blob=entry[2]
 assert git('diff','--name-only','--',name)=='' and git('diff','--cached','--name-only','--',name)==''
 with target.open('rb') as f: actual=digest(f)
 with subprocess.Popen(['git','-C',str(a.repo),'cat-file','blob',blob],stdout=subprocess.PIPE) as child:
  stored=digest(child.stdout);assert child.wait()==0
 assert actual==stored
 handles=subprocess.run(['/usr/sbin/lsof',str(target)],capture_output=True,text=True)
 assert handles.returncode==1 and not handles.stdout.strip() and not handles.stderr.strip()
 manifest.append(dict(path=name,bytes=target.stat().st_size,SHA256=actual,gitBlob=blob))
record=dict(head=head,repository=str(a.repo),files=manifest,beforeAvailableBytes=free(),verifiedUnix=time.time(),
 preservation='Exact complete file bytes retained in committed Git blobs; only reproducible worktree copies removed.',deleted=[])
a.record.write_text(json.dumps(record,indent=2)+'\n')
for item in manifest:
 target=a.repo/item['path']
 with target.open('rb') as f:assert digest(f)==item['SHA256']
 handles=subprocess.run(['/usr/sbin/lsof',str(target)],capture_output=True,text=True)
 assert handles.returncode==1 and not handles.stdout.strip() and not handles.stderr.strip()
 subprocess.run(['git','-C',str(a.repo),'update-index','--skip-worktree','--',item['path']],check=True)
 target.unlink();record['deleted'].append(item['path'])
 a.record.write_text(json.dumps(record,indent=2)+'\n')
assert git('rev-parse','HEAD')==head
for item in manifest:
 assert not (a.repo/item['path']).exists()
 assert git('cat-file','-s',item['gitBlob'])==str(item['bytes'])
record.update(afterAvailableBytes=free(),finishedUnix=time.time())
record['recoveredBytes']=record['afterAvailableBytes']-record['beforeAvailableBytes']
a.record.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
