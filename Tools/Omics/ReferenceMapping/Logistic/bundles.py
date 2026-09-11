"""Exact content-deduplicated tar bundles, with explicit reconstruction paths."""
import hashlib,io,json,tarfile,shutil,subprocess
from pathlib import Path,PurePosixPath
from prepare import sha

def pack(root,archive):
 files=sorted(p for p in root.rglob('*') if p.is_file());assert all(not p.is_symlink() for p in root.rglob('*'))
 members={str(p.relative_to(root)):dict(SHA256=sha(p),bytes=p.stat().st_size) for p in files};objects={}
 for p in files:objects.setdefault(members[str(p.relative_to(root))]['SHA256'],p)
 raw=(json.dumps(members,sort_keys=True,indent=2)+'\n').encode()
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  info=tarfile.TarInfo('members.json');info.size=len(raw);t.addfile(info,io.BytesIO(raw))
  for h,p in sorted(objects.items()):t.add(p,arcname='objects/'+h,recursive=False)
 verify(archive)
 return dict(members=members,bytes=archive.stat().st_size,SHA256=sha(archive),decodedBytes=sum(x['bytes'] for x in members.values()),uniqueBytes=sum(p.stat().st_size for p in objects.values()),objects=len(objects))
def verify(archive):
 with tarfile.open(archive,'r:gz') as t:
  items=t.getmembers();assert len({m.name for m in items})==len(items) and all(m.isfile() for m in items)
  members=json.load(t.extractfile('members.json'));objects={x['SHA256'] for x in members.values()}
  assert {m.name for m in items}=={'members.json'}|{'objects/'+h for h in objects}
  for path,v in members.items():
   parts=PurePosixPath(path);assert not parts.is_absolute() and '..' not in parts.parts and str(parts)==path
   assert len(v['SHA256'])==64 and all(c in '0123456789abcdef' for c in v['SHA256'])
   assert t.getmember('objects/'+v['SHA256']).size==v['bytes']
  # Follow archive order so gzip is decoded sequentially, not restarted for
  # each object in hash-set order. Every object is still checked in full.
  for item in items:
   if not item.name.startswith('objects/'):continue
   h=item.name.split('/')[1]
   digest=hashlib.sha256()
   with t.extractfile('objects/'+h) as f:
    for b in iter(lambda:f.read(1048576),b''):digest.update(b)
   assert digest.hexdigest()==h
 return members
def read_json(t,members,path):return json.load(t.extractfile('objects/'+members[path]['SHA256']))
def restore(archive,out):
 members=verify(archive);out.mkdir(parents=True,exist_ok=False)
 restored={}
 with tarfile.open(archive,'r:gz') as t:
  for name,v in members.items():
   p=out/name;p.parent.mkdir(parents=True,exist_ok=True)
   previous=restored.get(v['SHA256'])
   cloned=False
   if previous is not None:
    cloned=subprocess.run(['/bin/cp','-c',str(previous),str(p)],capture_output=True).returncode==0
   if not cloned:
    assert not p.exists() and shutil.disk_usage(out).free>v['bytes']+200000000
    with t.extractfile('objects/'+v['SHA256']) as f,p.open('xb') as g:
     for b in iter(lambda:f.read(1048576),b''):g.write(b)
   assert sha(p)==v['SHA256'] and p.stat().st_size==v['bytes']
   restored[v['SHA256']]=p
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('--restore',type=Path);a=p.parse_args()
 if a.restore:restore(a.archive,a.restore)
 else:print(json.dumps(dict(status='verified',members=len(verify(a.archive)))))
