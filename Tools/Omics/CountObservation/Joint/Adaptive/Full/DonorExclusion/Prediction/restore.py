"""Restore a retained evidence archive; optionally link qualified sparse caches."""
from pathlib import Path
import hashlib,json,os,shutil,sys,tarfile
source,dest=map(Path,sys.argv[1:3]);dest.mkdir(exist_ok=False);index=json.loads((source/'files.json').read_text());receipt=json.loads((source/'retention.json').read_text());assert hashlib.sha256((source/'files.json').read_bytes()).hexdigest()==receipt['filesIndexSHA256'];objects=dest/'objects';objects.mkdir();expected={v['sha256']:v['bytes'] for v in index['files'].values()}
class Reader:
 def __init__(self):self.i=0;self.f=None;self.hash=hashlib.sha256()
 def read(self,n=-1):
  assert n>=0;out=b''
  while len(out)<n and self.i<len(receipt['parts']):
   part=receipt['parts'][self.i]
   if self.f is None:
    p=source/part['path'];assert p.name==part['path'] and p.stat().st_size==part['bytes'];assert hashlib.sha256(p.read_bytes()).hexdigest()==part['sha256'];self.f=p.open('rb')
   b=self.f.read(n-len(out));out+=b;self.hash.update(b)
   if not b:self.f.close();self.f=None;self.i+=1
  return out
reader=Reader();seen=set()
with tarfile.open(fileobj=reader,mode='r|') as tar:
 for member in tar:
  digest=member.name.removeprefix('objects/');assert member.name=='objects/'+digest and digest in expected and digest not in seen and member.isfile() and member.size==expected[digest]
  h=hashlib.sha256()
  with (objects/digest).open('xb') as dst:
   src=tar.extractfile(member)
   while b:=src.read(1<<20):h.update(b);dst.write(b)
  assert h.hexdigest()==digest;seen.add(digest)
while reader.read(1<<20):pass
assert reader.hash.hexdigest()==receipt['archiveSHA256'] and seen==set(expected)
for name,e in index['files'].items():
 p=dest/name;assert not Path(name).is_absolute() and '..' not in Path(name).parts;p.parent.mkdir(parents=True,exist_ok=True)
 if p.name in ['verification.json','scores-verification.json']:shutil.copyfile(objects/e['sha256'],p)
 else:os.link(objects/e['sha256'],p)
if len(sys.argv)>3:
 full=Path(sys.argv[3])
 for tag in index['tags']:
  origin=tag.split('-')[0];p=full/(origin+'-cells.h5');import gzip
  m=json.loads(gzip.decompress((dest/'fits'/tag/(origin+'-manifest.json.gz')).read_bytes()));h=hashlib.sha256()
  with p.open('rb') as f:
   while b:=f.read(1<<20):h.update(b)
  assert h.hexdigest()==m['cacheSHA256'];os.link(p,dest/'fits'/tag/p.name)
print({'status':'restored-all-logical-files','files':len(index['files']),'objects':len(seen),'tags':index['tags']})
