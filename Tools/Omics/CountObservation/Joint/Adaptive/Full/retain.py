"""Stream completed full-gene evidence into bounded, content-addressed parts.

Each part is at most 48 MiB. Packing and verification never read the whole
archive or the sparse count cache into memory. No source files are removed.
"""
from pathlib import Path, PurePosixPath
import gzip,hashlib,json,os,shutil,sys,tarfile
CHUNK=1<<20
PART=48<<20

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while b:=f.read(CHUNK):h.update(b)
 return h.hexdigest()

class Parts:
 def __init__(self,root):self.root=root;self.current=None;self.parts=[];self.total=0;self.joined=hashlib.sha256()
 def write(self,data):
  size=len(data);view=memoryview(data)
  while view:
   if self.current is None:
    self.name=f'objects.tar.gz.part{len(self.parts):03d}';self.current=(self.root/(self.name+'.partial')).open('xb');self.used=0;self.digest=hashlib.sha256()
   chunk=view[:PART-self.used];self.current.write(chunk);self.digest.update(chunk);self.joined.update(chunk);self.used+=len(chunk);self.total+=len(chunk);view=view[len(chunk):]
   if self.used==PART:self.finish_part()
  return size
 def finish_part(self):
  if self.current is None:return
  self.current.close();self.current=None;(self.root/(self.name+'.partial')).rename(self.root/self.name)
  self.parts.append({'name':self.name,'bytes':self.used,'SHA256':self.digest.hexdigest()})
 def flush(self):
  if self.current:self.current.flush()
 def finish(self):
  self.finish_part();return {'bytes':self.total,'SHA256':self.joined.hexdigest(),'parts':self.parts}

class Reader:
 def __init__(self,paths):self.paths=iter(paths);self.current=None
 def read(self,n):
  b=bytearray()
  while len(b)<n:
   if self.current is None:
    try:self.current=next(self.paths).open('rb')
    except StopIteration:break
   piece=self.current.read(n-len(b))
   if piece:b.extend(piece)
   else:self.current.close();self.current=None
  return bytes(b)
 def close(self):
  if self.current:self.current.close()

def unpack(evidence,destination=None):
 m=json.loads((evidence/'manifest.json').read_bytes());joined=hashlib.sha256();total=0;paths=[]
 for p in m['archive']['parts']:
  assert PurePosixPath(p['name']).name==p['name'];path=evidence/p['name'];h=hashlib.sha256();n=0
  with path.open('rb') as f:
   while b:=f.read(CHUNK):h.update(b);joined.update(b);n+=len(b)
  assert n==p['bytes'] and h.hexdigest()==p['SHA256'];total+=n;paths.append(path)
 assert total==m['archive']['bytes'] and joined.hexdigest()==m['archive']['SHA256']
 aliases={}
 for n,v in m['files'].items():
  p=PurePosixPath(n);assert not p.is_absolute() and '..' not in p.parts
  assert v['SHA256'] in m['objects'] and v['bytes']==m['objects'][v['SHA256']]['bytes'];aliases.setdefault(v['SHA256'],[]).append((n,v))
 assert set(aliases)==set(m['objects'])
 if destination:destination.mkdir(parents=True,exist_ok=False)
 reader=Reader(paths);seen=set()
 try:
  with tarfile.open(fileobj=reader,mode='r|gz') as tar:
   for member in tar:
    assert member.isfile() and member.name.startswith('objects/');key=member.name[8:];assert key in m['objects'] and key not in seen;seen.add(key)
    assert member.size==m['objects'][key]['bytes'];h=hashlib.sha256();n=0;out=None
    if destination:
     first=destination/aliases[key][0][0];first.parent.mkdir(parents=True,exist_ok=True);out=first.open('xb')
    try:
     with tar.extractfile(member) as f:
      while b:=f.read(CHUNK):h.update(b);n+=len(b);out.write(b) if out else None
    finally:
     if out:out.close()
    assert n==member.size and h.hexdigest()==key
    if destination:
     first.chmod(aliases[key][0][1]['mode'])
     for name,v in aliases[key][1:]:
      p=destination/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(first,p);p.chmod(v['mode'])
 finally:reader.close()
 assert seen==set(m['objects']);print(json.dumps({'verifiedLogicalFiles':len(m['files']),'verifiedObjects':len(seen),'restored':destination is not None}))

def pack(study,repo,out,origins):
 assert origins and all(o in ['Kang','HIRISA'] for o in origins)
 out.mkdir(parents=True,exist_ok=False);files={};objects={};signatures={}
 def add(name,p):
  assert p.is_file() and not p.is_symlink();s=p.stat();h=sha(p);assert (s.st_ino,s.st_size,s.st_mtime_ns)==(p.stat().st_ino,p.stat().st_size,p.stat().st_mtime_ns)
  files[name]={'SHA256':h,'bytes':s.st_size,'mode':0o755 if os.access(p,os.X_OK) else 0o644};objects.setdefault(h,p);signatures[p]=(s.st_ino,s.st_size,s.st_mtime_ns)
 for origin in origins:
  folder=study/(origin+'-gap-repaired');v=json.loads((folder/'verification-state.json').read_text());w=json.loads((folder/'verification-watch.json').read_text())
  assert w['status']=='completed-all-workers-and-verification' and not v['errors'];assert v['genesVerified']==json.loads((study/(origin+'-prepare-state.json')).read_text())['genes']
  for name in [origin,origin+'-gap-repaired']:
   for p in sorted((study/name).iterdir()):
    assert '.partial' not in p.name
    if p.is_file() and not p.name.endswith('.lock'):add('study/'+name+'/'+p.name,p)
  for suffix in ['-cells.h5','-manifest.json.gz','-prepare-state.json','-prepare.log']:
   add('study/'+origin+suffix,study/(origin+suffix))
 for p in sorted(study.iterdir()):
  if p.is_file() and (p.name in ['protocol.json','results.json','runtime-freeze.json','full-counts','baseline-run.py','build.log'] or p.name.endswith('recipe-freeze.json') or p.name=='full-recipes-freeze.json'):
   add('study/'+p.name,p)
 for p in sorted((study/'gap-repair').iterdir()):
  if p.is_file():add('study/gap-repair/'+p.name,p)
 calibration=os.environ.get('NUMIVIVO_DONOR_EXCLUSION_STUDY')
 if calibration:
  cal=Path(calibration);state=json.loads((cal/'state.json').read_text());checked=json.loads((cal/'verification.json').read_text())
  assert state['status']=='completed-all-training-only-folds' and checked['status']=='passed-all-training-only-calibrations'
  for p in sorted(cal.rglob('*')):
   if p.is_file() and '__pycache__' not in p.parts and 'publication' not in p.relative_to(cal).parts:
    assert '.partial' not in p.name;add('donor-exclusion/study/'+str(p.relative_to(cal)),p)
  for p in sorted((repo/'Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion').iterdir()):
   if p.is_file():add('donor-exclusion/recipes/'+p.name,p)
  for origin in ['Kang','HIRISA']:
   add('study/'+origin+'-manifest.json.gz',study/(origin+'-manifest.json.gz'))
  add('source/VivoCountObservationCalibration.swift',repo/'Sources/NumiVivoKit/Omics/VivoCountObservationCalibration.swift')
 for p in sorted((repo/'Tools/Omics/CountObservation/Joint/Adaptive/Full').iterdir()):
  if p.is_file():add('recipes/'+p.name,p)
 for name in ['VivoCountObservation.swift','VivoCountRateLikelihood.swift','VivoJointCountResponse.swift','VivoAdaptiveJointCountResponse.swift']:
  add('source/'+name,repo/'Sources/NumiVivoKit/Omics'/name)
 add('source/SingleCellJointCountTests.swift',repo/'Tests/NumiVivoIntegrationTests/SingleCellJointCountTests.swift')
 writer=Parts(out)
 with gzip.GzipFile(fileobj=writer,filename='',mode='wb',mtime=0) as gz,tarfile.open(fileobj=gz,mode='w|') as tar:
  for h,p in sorted(objects.items()):
   s=p.stat();assert (s.st_ino,s.st_size,s.st_mtime_ns)==signatures[p];info=tarfile.TarInfo('objects/'+h);info.size=s.st_size;info.mode=0o644;info.mtime=0
   with p.open('rb') as f:tar.addfile(info,f)
   assert sha(p)==h
 archive=writer.finish();m={'schemaVersion':1,'format':'sha256-objects-with-logical-file-aliases','origins':origins,'files':files,'objects':{h:{'bytes':p.stat().st_size} for h,p in objects.items()},'archive':archive,'scope':'Completed origin full-gene original-cell cache, baseline and repaired native fits including failures, independent checks, frozen executable and source, protocols and recipes. Training fitting only; no new outcome validation.'}
 if 'HIRISA' in origins:
  prior=repo/'Tools/Omics/CountObservation/Joint/Adaptive/Full/evidence/2026-09-11-kang/manifest.json'
  prior_manifest=json.loads(prior.read_text())
  m['relatedPublishedEvidence']={'Kang':{'commit':'217614c515e0666cf795ffb9d94b31a49789129e','path':str(prior.relative_to(repo).parent),'archiveSHA256':prior_manifest['archive']['SHA256'],'note':'Kang complete raw fits remain in this separate archive; both full-cohort result summaries may be reported together.'}}
 (out/'manifest.json').write_text(json.dumps(m,indent=2,sort_keys=True)+'\n');(out/'results.json').write_bytes((study/'results.json').read_bytes());print(json.dumps({'files':len(files),'objects':len(objects),'archiveBytes':archive['bytes'],'archiveSHA256':archive['SHA256']}))

if __name__=='__main__':
 if sys.argv[1]=='verify':unpack(Path(sys.argv[2]))
 elif sys.argv[1]=='restore':unpack(Path(sys.argv[2]),Path(sys.argv[3]))
 elif sys.argv[1]=='pack':pack(*map(Path,sys.argv[2:5]),sys.argv[5:])
 else:raise ValueError('Use pack STUDY REPO OUTPUT ORIGIN..., verify EVIDENCE, or restore EVIDENCE DESTINATION')
