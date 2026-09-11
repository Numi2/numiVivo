"""Retain streamed-report qualification, reusing exact objects in the prior archive."""
import argparse,hashlib,io,json,os,shutil,subprocess,tarfile
from pathlib import Path,PurePosixPath
from run_parse import sha,save

def git_bytes(repo,commit,path):return subprocess.check_output(['git','-C',str(repo),'cat-file','blob',commit+':'+path])
def validate(archive):
 with tarfile.open(archive,'r:gz') as t:
  items=t.getmembers();assert all(i.isfile() for i in items) and len({i.name for i in items})==len(items)
  members=json.load(t.extractfile('members.json'));dependencies=json.load(t.extractfile('dependencies.json'))
  for name,v in members.items():
   p=PurePosixPath(name);assert not p.is_absolute() and '..' not in p.parts and str(p)==name
   assert v['archive'] in ('self','baseline') and v['bytes']>=0 and 0<=v['mode']<=511
   assert len(v['SHA256'])==64 and all(c in '0123456789abcdef' for c in v['SHA256'])
  own={v['SHA256'] for v in members.values() if v['archive']=='self'}
  assert {i.name for i in items}=={'members.json','dependencies.json'}|{'objects/'+h for h in own}
  for v in members.values():
   if v['archive']=='self':assert t.getmember('objects/'+v['SHA256']).size==v['bytes']
  for h in own:
   digest=hashlib.sha256()
   with t.extractfile('objects/'+h) as f:
    for b in iter(lambda:f.read(1048576),b''):digest.update(b)
   assert digest.hexdigest()==h
 return members,dependencies

def restore(a):
 members,dependency=validate(a.archive);assert not a.out.exists();a.out.mkdir(parents=True);counts=dict(sourceAPFSCopies=0,withinDestinationAPFSCopies=0,directExtractions=0,verifiedFileCopies=0);done=set();byhash={}
 def place(t,item,origin):
  digest=item.name.removeprefix('objects/');targets=[(n,v) for n,v in members.items() if v['archive']==origin and v['SHA256']==digest]
  for n,v in targets:
   assert item.isfile() and item.size==v['bytes'];p=a.out/n;p.parent.mkdir(parents=True,exist_ok=True);previous=byhash.get(digest);candidate=previous or (a.study/n if a.study else None);copied=False
   if candidate and candidate.is_file() and not candidate.is_symlink() and sha(candidate)==digest:
    copied=subprocess.run(['/bin/cp','-c',str(candidate),str(p)],capture_output=True).returncode==0
    if copied:assert p.stat().st_ino!=candidate.stat().st_ino;counts['withinDestinationAPFSCopies' if previous else 'sourceAPFSCopies']+=1
   if not copied:
    assert not p.exists() and shutil.disk_usage(a.out).free>v['bytes']+16000000
    with (previous.open('rb') if previous else t.extractfile(item)) as f,p.open('xb') as g:shutil.copyfileobj(f,g,1048576)
    counts['verifiedFileCopies' if previous else 'directExtractions']+=1
   assert sha(p)==digest;p.chmod(v['mode']);done.add(n);byhash[digest]=p
 with tarfile.open(a.archive,'r:gz') as t:
  for item in t.getmembers():
   if item.name.startswith('objects/'):place(t,item,'self')
 process=subprocess.Popen(['git','-C',str(a.repo),'cat-file','blob',dependency['commit']+':'+dependency['path']],stdout=subprocess.PIPE)
 class Reader:
  def __init__(self):self.h=hashlib.sha256()
  def read(self,n=-1):
   b=process.stdout.read(n);self.h.update(b);return b
 reader=Reader()
 with tarfile.open(fileobj=reader,mode='r|gz') as t:
  for item in t:
   if item.name.startswith('objects/'):place(t,item,'baseline')
 while reader.read(1048576):pass
 assert process.wait()==0 and reader.h.hexdigest()==dependency['SHA256'] and done==set(members)
 result=dict(status='restored-all-qualified-bytes',members=len(members),counts=counts,archiveSHA256=sha(a.archive),baselineArchiveSHA256=dependency['SHA256'],allHashesExact=True,nativeReexecutedDuringRestoration=False,predictionFitOrScoring=False)
 save(a.out/'restoration.json',result);print(json.dumps(result))

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=['pack','verify','restore']);p.add_argument('--study',type=Path);p.add_argument('--repo',type=Path);p.add_argument('--out',type=Path);p.add_argument('--archive',type=Path);a=p.parse_args()
 if a.action=='verify':m,d=validate(a.archive);print(json.dumps(dict(status='verified-self-archive',members=len(m),requiredBaseline=d)));return
 if a.action=='restore':restore(a);return
 s=a.study.resolve();r=a.repo.resolve();result=json.loads((s/'pipeline.json').read_text());assert result['status']=='passed-complete-native-reconstruction-and-comparison'
 freeze=json.loads((s/'source-freeze.json').read_text())
 for name,h in freeze['files'].items():assert sha(Path(name))==h,name
 protocol=json.loads((s/'protocol.json').read_text());baseline=json.loads(git_bytes(r,protocol['baselineCommit'],str(Path(protocol['baselineArchivePath']).parent/'contents.json')))
 assert baseline['archiveSHA256']==protocol['baselineArchiveSHA256'];known={v['SHA256']:v['bytes'] for v in baseline['members'].values()}
 dependency=dict(commit=protocol['baselineCommit'],path=protocol['baselineArchivePath'],SHA256=protocol['baselineArchiveSHA256'])
 files={}
 for p in (s/'native-file').rglob('*'):
  if p.is_file():files[str(p.relative_to(s))]=p
 for name in ['protocol.json','plan.json','source-freeze.json','pipeline.json','publish.json','verify.json','results.json','profile-provenance.json','Profile.swift','profile','profile.log','profile.stderr','build.log','tests.log','pipeline.log','publish.stdout','publish.stderr','verify.stdout','verify.stderr','retention-preflight.json','retention-preflight.log','prior-archive-dehydration.json']:
  p=s/name;assert p.is_file(),name;files[name]=p
 for n in ['build/numivivo-omics','build/sources.sha256']:files[n]=s/n
 for name in freeze['files']:
  p=Path(name)
  if p.is_relative_to(r):files['implementation/'+str(p.relative_to(r))]=p
 for p in Path(__file__).parent.iterdir():
  if p.is_file() and p.suffix in ('.py','.swift'):files['retention-recipe/'+p.name]=p
 assert all(not p.is_symlink() for p in files.values());members={};objects={}
 for n,p in sorted(files.items()):
  h=sha(p);size=p.stat().st_size;origin='baseline' if known.get(h)==size else 'self';members[n]=dict(SHA256=h,bytes=size,mode=p.stat().st_mode&511,archive=origin)
  if origin=='self':objects.setdefault(h,p)
 assert not a.out.exists();temp=a.out.with_name(a.out.name+'.incomplete');assert not temp.exists();temp.mkdir(parents=True);archive=temp/'results.tar.gz'
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  for n,v in [('members.json',members),('dependencies.json',dependency)]:
   raw=(json.dumps(v,sort_keys=True,indent=2)+'\n').encode();item=tarfile.TarInfo(n);item.size=len(raw);t.addfile(item,io.BytesIO(raw))
  for h,p in sorted(objects.items()):
   item=tarfile.TarInfo('objects/'+h);item.size=p.stat().st_size
   with p.open('rb') as f:t.addfile(item,f)
 assert validate(archive)==(members,dependency)
 for n,p in files.items():assert sha(p)==members[n]['SHA256']
 save(temp/'contents.json',dict(status='verified-self-archive',archiveSHA256=sha(archive),archiveBytes=archive.stat().st_size,members=members,requiredBaseline=dependency,scope='All native output bytes, new executable/source, phase profile, tests, whole-cohort metrics and exact-byte comparison. Identical objects are retained in the exact earlier archive; no duplicate report or count-source payload. Both archives are needed to restore. Parse derivatives CC BY-NC 4.0.'))
 for n in ['pipeline.json','protocol.json','source-freeze.json']:shutil.copyfile(s/n,temp/n)
 records=[]
 for p in sorted(temp.iterdir()):
  records.append(dict(sourcePath=p.name,sourceBytes=p.stat().st_size,sourceSHA256=sha(p),storedPath=p.name,storedBytes=p.stat().st_size,storedSHA256=sha(p),gzipEncoded=False,bundledSourceFiles=False))
 save(temp/'manifest.json',dict(schemaVersion=1,records=records));os.rename(temp,a.out);print(json.dumps(dict(status='retained-complete-streamed-report-evidence',members=len(members),archiveBytes=(a.out/'results.tar.gz').stat().st_size)))
if __name__=='__main__':main()
