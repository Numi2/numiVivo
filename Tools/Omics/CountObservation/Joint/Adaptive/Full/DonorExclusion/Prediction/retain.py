"""Retain completed fold evidence with deduplicated objects and bounded parts."""
from pathlib import Path
import hashlib,io,json,sys,tarfile
study,pred,dest=map(Path,sys.argv[1:4]);tags=sys.argv[4:];assert tags;dest.mkdir(parents=True,exist_ok=False)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()
files={}
def add(path,name):assert path.is_file() and name not in files;files[name]={'sha256':sha(path),'bytes':path.stat().st_size,'source':str(path)}
for p in Path(__file__).parent.iterdir():
 if p.is_file() and p.suffix in ['.py','.swift','.sh','.md']:add(p,'recipes/'+p.name)
for name in ['runtime-freeze.json','recipe-freeze.json','scoring-plan.json','scoring-initial-failure.json','Run.swift','run.py','verify.py','reference.py','score.py','verify_scores.py','check_contract.py','advance.py','predict','build.log']:
 add(pred/name,'prediction/'+name)
for p in (pred/'contract').iterdir():
 if p.is_file():add(p,'prediction/contract/'+p.name)
for tag in tags:
 assert json.loads((pred/tag/'scores-verification.json').read_text())['status']=='passed-independent-sparse-scoring'
 origin=tag.split('-')[0]
 for p in (pred/tag).iterdir():
  if p.is_file() and p.suffix not in ['.lock','.tmp']:add(p,'prediction/'+tag+'/'+p.name)
 for p in (study/tag/origin).iterdir():
  if p.is_file() and p.suffix not in ['.lock','.tmp']:add(p,'fits/'+tag+'/'+origin+'/'+p.name)
 for name in [origin+'-manifest.json.gz',origin+'-prepare-state.json','runtime-freeze.json']:add(study/tag/name,'fits/'+tag+'/'+name)
add(study/'protocol.json','fits/protocol.json')
objects={}
for entry in files.values():objects.setdefault(entry['sha256'],entry)
index={'schemaVersion':1,'tags':tags,'qualification':'Completed reused development donor folds only; other folds remain separate active work. No calibrated interval or independent biological validation claim.','cacheRestoration':'Restore Kang-cells.h5 or HIRISA-cells.h5 from the published full-cohort evidence archives, verify cacheSHA256 against each fold manifest, then hard-link into fits/TAG/. No source cache is duplicated here.','files':{k:{a:v[a] for a in ['sha256','bytes']} for k,v in files.items()}}
(dest/'files.json').write_text(json.dumps(index,indent=2)+'\n');parts=[]
class Writer:
 def __init__(self):self.file=None;self.n=0;self.total=0;self.hash=hashlib.sha256();self.path=None
 def write(self,b):
  self.hash.update(b);self.total+=len(b)
  while b:
   if self.file is None:self.path=dest/f'evidence.tar.part-{len(parts):03d}';self.file=self.path.open('xb');self.n=0
   take=min(len(b),48*1024*1024-self.n);self.file.write(b[:take]);self.n+=take;b=b[take:]
   if self.n==48*1024*1024:self.finish()
  return None
 def finish(self):
  if self.file:self.file.close();parts.append({'path':self.path.name,'bytes':self.n,'sha256':sha(self.path)});self.file=None
 def flush(self):
  if self.file:self.file.flush()
w=Writer()
with tarfile.open(fileobj=w,mode='w|') as tar:
 for digest,e in sorted(objects.items()):
  p=Path(e['source']);assert sha(p)==digest;info=tarfile.TarInfo('objects/'+digest);info.size=e['bytes'];info.mode=0o644;info.mtime=0
  with p.open('rb') as f:tar.addfile(info,f)
w.finish()
# Stream all parts back through tar and verify every logical object.
class Reader:
 def __init__(self):self.i=0;self.file=None
 def read(self,n=-1):
  assert n>=0;out=b''
  while len(out)<n and self.i<len(parts):
   if self.file is None:self.file=(dest/parts[self.i]['path']).open('rb')
   b=self.file.read(n-len(out));out+=b
   if not b:self.file.close();self.file=None;self.i+=1
  return out
seen=set()
with tarfile.open(fileobj=Reader(),mode='r|') as tar:
 for member in tar:
  digest=member.name.removeprefix('objects/');assert digest in objects and digest not in seen;f=tar.extractfile(member);h=hashlib.sha256();size=0
  while b:=f.read(1<<20):h.update(b);size+=len(b)
  assert h.hexdigest()==digest and size==objects[digest]['bytes'];seen.add(digest)
assert seen==set(objects)
r={'status':'retained-and-all-objects-readback-verified','parts':parts,'archiveBytes':w.total,'archiveSHA256':w.hash.hexdigest(),'logicalFiles':len(files),'uniqueObjects':len(objects),'filesIndexSHA256':sha(dest/'files.json'),'tags':tags};(dest/'retention.json').write_text(json.dumps(r,indent=2)+'\n');print(r)
