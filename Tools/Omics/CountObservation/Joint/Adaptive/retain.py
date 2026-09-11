"""Content-addressed, lossless retention with aliases for identical complete replays."""
from pathlib import Path,PurePosixPath
import gzip,hashlib,io,json,os,sys,tarfile

def digest(b):return hashlib.sha256(b).hexdigest()
class PartReader:
 def __init__(self,paths):self.paths=iter(paths);self.current=None
 def read(self,n):
  data=bytearray()
  while len(data)<n:
   if self.current is None:
    try:self.current=next(self.paths).open('rb')
    except StopIteration:break
   chunk=self.current.read(n-len(data))
   if chunk:data.extend(chunk)
   else:self.current.close();self.current=None
  return bytes(data)
 def close(self):
  if self.current:self.current.close()

def restore(evidence,destination):
 manifest=json.loads((evidence/'manifest.json').read_bytes());paths=[];joined=hashlib.sha256()
 for part in manifest['archive']['parts']:
  name=part['name'];assert PurePosixPath(name).name==name
  p=evidence/name;raw=p.read_bytes();assert len(raw)==part['bytes'] and digest(raw)==part['SHA256'];joined.update(raw);paths.append(p)
 assert joined.hexdigest()==manifest['archive']['SHA256']
 destination.mkdir(parents=True,exist_ok=False);aliases={};seen=set()
 for name,entry in manifest['files'].items():
  path=PurePosixPath(name);assert not path.is_absolute() and '..' not in path.parts
  aliases.setdefault(entry['SHA256'],[]).append((name,entry))
 reader=PartReader(paths)
 try:
  with tarfile.open(fileobj=reader,mode='r|gz') as tar:
   for member in tar:
    assert member.isfile() and member.name.startswith('objects/')
    h=member.name[len('objects/'):];assert h in manifest['objects'] and h not in seen;seen.add(h)
    raw=tar.extractfile(member).read();assert len(raw)==manifest['objects'][h]['bytes'] and digest(raw)==h
    for name,entry in aliases[h]:
     assert len(raw)==entry['bytes'];output=destination/name;output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(raw);output.chmod(entry['mode'])
 finally:reader.close()
 assert seen==set(manifest['objects'])
 print(json.dumps({'restoredLogicalFiles':len(manifest['files']),'verifiedUniqueObjects':len(seen)}))

if sys.argv[1]=='restore':restore(Path(sys.argv[2]),Path(sys.argv[3]));sys.exit()
study,repo,out=map(Path,sys.argv[1:4]);out.mkdir(parents=True,exist_ok=False)
files={};objects={}
def retain(name,p):
 raw=p.read_bytes();h=digest(raw);objects[h]=raw
 files[name]={'SHA256':h,'bytes':len(raw),'mode':0o755 if os.access(p,os.X_OK) else 0o644}
for p in sorted(study.rglob('*')):
 if not p.is_file() or '__pycache__' in p.parts:continue
 rel=p.relative_to(study)
 if rel.parts[0]=='publication' or rel.parts[0].startswith('restored'):continue
 retain('study/'+str(rel),p)
for p in sorted((repo/'Tools/Omics/CountObservation/Joint/Adaptive').iterdir()):
 if p.is_file():retain('recipes/'+p.name,p)
for name in ['VivoCountObservation.swift','VivoCountRateLikelihood.swift','VivoJointCountResponse.swift','VivoAdaptiveJointCountResponse.swift']:
 retain('source/'+name,repo/'Sources/NumiVivoKit/Omics'/name)
retain('source/SingleCellJointCountTests.swift',repo/'Tests/NumiVivoIntegrationTests/SingleCellJointCountTests.swift')
parent_manifest=repo/'Tools/Omics/CountObservation/Joint/evidence/2026-09-11/manifest.json'
retain('reference/parent-manifest.json',parent_manifest)
parent=Path('/Users/home/numivivo-joint-counts-20260911')
for origin in ['Kang','HIRISA']:
 retain('study/parent/'+origin+'-native.jsonl.gz',parent/(origin+'-native.jsonl.gz'))
external={'parentFiniteGridArchive':{'commit':'5cc9bde1d3c330657030426c85750720eeac2743','path':'Tools/Omics/CountObservation/Joint/evidence/2026-09-11','SHA256':json.loads(parent_manifest.read_bytes())['archive']['SHA256'],'purpose':'Previous finite-grid qualification, source lineage and retained failures. Necessary parent comparison outputs are included in this archive.'}}
archive=out/'objects.tar.gz'
with archive.open('wb') as f,gzip.GzipFile(fileobj=f,filename='',mode='wb',mtime=0) as gz,tarfile.open(fileobj=gz,mode='w') as tar:
 for h,raw in sorted(objects.items()):
  info=tarfile.TarInfo('objects/'+h);info.size=len(raw);info.mode=0o644;info.mtime=0;tar.addfile(info,io.BytesIO(raw))
parts=[];archive_bytes=archive.stat().st_size;archive_hash=digest(archive.read_bytes())
with archive.open('rb') as f:
 while raw:=f.read(48*1024*1024):
  p=out/('objects.tar.gz.part%03d'%len(parts));p.write_bytes(raw);parts.append({'name':p.name,'bytes':len(raw),'SHA256':digest(raw)})
assert sum(p['bytes'] for p in parts)==archive_bytes
joined=hashlib.sha256()
for part in parts:joined.update((out/part['name']).read_bytes())
assert joined.hexdigest()==archive_hash
archive.unlink()
manifest={'schemaVersion':1,'format':'sha256-objects-with-logical-file-aliases','files':files,'objects':{h:{'bytes':len(b)} for h,b in objects.items()},'externalSources':external,'archive':{'bytes':archive_bytes,'SHA256':archive_hash,'parts':parts},'scope':'Complete adaptive initial and final fits, support partitions, query predictions, independent per-cell references, previous grid65 comparison outputs, exact binaries and source snapshots, protocols, logs and retained iteration-limit failure. Original source lineage is bound to the published parent archive.'}
(out/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n');(out/'results.json').write_bytes((study/'results.json').read_bytes())
print(json.dumps({'logicalFiles':len(files),'uniqueObjects':len(objects),'archive':manifest['archive']}))
