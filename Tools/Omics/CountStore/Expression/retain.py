"""Retain complete Parse expression evidence; restore verified separate-inode files."""
import argparse,gzip,hashlib,io,json,os,shutil,subprocess,tarfile,time
from pathlib import Path,PurePosixPath
from report_delta import sha,chunks

def write(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
def verify(p):
 with tarfile.open(p,'r:gz') as t:
  items=t.getmembers();assert len({i.name for i in items})==len(items) and all(i.isfile() for i in items)
  m=json.load(t.extractfile('members.json'));assert {i.name for i in items}=={'members.json'}|{'objects/'+v['SHA256'] for v in m.values()}
  for name,v in m.items():
   n=PurePosixPath(name);h=v['SHA256'];assert str(n)==name and not n.is_absolute() and '..' not in n.parts
   assert len(h)==64 and all(c in '0123456789abcdef' for c in h) and 0<=v['mode']<=511
   assert t.getmember('objects/'+h).size==v['bytes']
  for i in items:
   if i.name=='members.json':continue
   h=hashlib.sha256()
   with t.extractfile(i) as f:
    for b in iter(lambda:f.read(1048576),b''):h.update(b)
   assert h.hexdigest()==i.name[8:]
 return m

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=['pack','verify','restore']);p.add_argument('--study',type=Path);p.add_argument('--repo',type=Path);p.add_argument('--archive',type=Path);p.add_argument('--out',type=Path);a=p.parse_args()
 if a.action=='verify':print(json.dumps(dict(status='verified',members=len(verify(a.archive)))));return
 if a.action=='restore':
  members=verify(a.archive);assert not a.out.exists();a.out.mkdir(parents=True);clones=0;extractions=0
  with tarfile.open(a.archive,'r:gz') as t:
   for n,v in members.items():
    target=a.out/n;target.parent.mkdir(parents=True,exist_ok=True);source=a.study/n if a.study else None;cloned=False
    if source and source.is_file() and not source.is_symlink() and sha(source)==v['SHA256']:
     cloned=subprocess.run(['/bin/cp','-c',str(source),str(target)],capture_output=True).returncode==0
     if cloned:assert source.stat().st_ino!=target.stat().st_ino;clones+=1
    if not cloned:
     assert not target.exists() and shutil.disk_usage(a.out).free>v['bytes']+16000000
     with t.extractfile('objects/'+v['SHA256']) as f,target.open('xb') as g:shutil.copyfileobj(f,g,1048576)
     extractions+=1
    assert sha(target)==v['SHA256'];target.chmod(v['mode'])
  descriptor=json.loads((a.out/'baseline-delta/descriptor.json').read_text());h=hashlib.sha256();size=0
  for b in chunks(a.out/'native-file/report.json',a.out/'baseline-delta'):h.update(b);size+=len(b)
  assert h.hexdigest()==descriptor['nativeReportSHA256'] and size==descriptor['nativeReportBytes']
  result=dict(status='restored-archive-and-exact-baseline',members=len(members),archiveSHA256=sha(a.archive),sourceAPFSClones=clones,directExtractions=extractions,baselineNativeSHA256=h.hexdigest(),nativeReplay='pending',predictionFitOrScoring=False)
  write(a.out/'restoration.json',result);print(json.dumps(result));return
 s=a.study.resolve();r=a.repo.resolve()
 assert json.loads((s/'native-pipeline.json').read_text())['status']=='passed-full-parse-native-aggregate-handoff'
 assert json.loads((s/'reference-pipeline.json').read_text())['status']=='completed-six-full-parse-reference-comparisons'
 assert json.loads((s/'baseline-delta/verification.json').read_text())['originalBytesComparedDirectly']
 for name in ['source-freeze.json','reference-freeze.json']:
  freeze=json.loads((s/name).read_text())
  for path,digest in freeze['files'].items():assert sha(Path(path))==digest,path
 assert not a.out.exists();temporary=a.out.with_name(a.out.name+'.incomplete');assert not temporary.exists();temporary.mkdir(parents=True)
 files={}
 for directory in ['native-file','baseline-delta','reference','reference-inputs','recipe']:
  for f in (s/directory).rglob('*'):
   if f.is_file() and '__pycache__' not in f.parts:files[str(f.relative_to(s))]=f
 for f in s.iterdir():
  if f.is_file() and (f.suffix in {'.json','.log','.tsv','.stdout','.stderr'} or f.name in ['native-reference']):files[f.name]=f
 for n in ['build/numivivo-omics','build/sources.sha256']:files[n]=s/n
 freeze=json.loads((s/'source-freeze.json').read_text())
 for path in freeze['files']:
  f=Path(path)
  if f.is_relative_to(r):files['implementation/'+str(f.relative_to(r))]=f
 for f in Path(__file__).parent.iterdir():
  if f.is_file():files['retention-recipe/'+f.name]=f
 assert all(not f.is_symlink() for f in files.values())
 members={n:dict(SHA256=sha(f),bytes=f.stat().st_size,mode=f.stat().st_mode&511) for n,f in sorted(files.items())};objects={}
 for n,f in files.items():objects.setdefault(members[n]['SHA256'],f)
 archive=temporary/'results.tar.gz';raw=(json.dumps(members,sort_keys=True,indent=2)+'\n').encode()
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  info=tarfile.TarInfo('members.json');info.size=len(raw);t.addfile(info,io.BytesIO(raw))
  for h,f in sorted(objects.items()):
   info=tarfile.TarInfo('objects/'+h);info.size=f.stat().st_size
   with f.open('rb') as source:t.addfile(info,source)
 assert verify(archive)==members
 assert all(sha(f)==members[n]['SHA256'] for n,f in files.items())
 write(temporary/'contents.json',dict(status='verified',members=members,archiveSHA256=sha(archive),archiveBytes=archive.stat().st_size,logicalBytes=sum(v['bytes'] for v in members.values()),scope='Complete native analysis bundle and owned count source, exact baseline JSON via lossless membership delta, both executables, frozen native source, all six R workflows and all-feature comparisons. Original baseline gzip transport encoding is not reproduced. Source counts were qualified separately; no prediction fitting or scoring. Parse derivatives CC BY-NC 4.0.'))
 for n in ['native-comparison.json','reference-comparison.json','protocol.json','source-freeze.json']:shutil.copyfile(s/n,temporary/n)
 records=[]
 for f in sorted(temporary.iterdir()):
  records.append(dict(sourcePath=f.name,sourceBytes=f.stat().st_size,sourceSHA256=sha(f),storedPath=f.name,storedBytes=f.stat().st_size,storedSHA256=sha(f),gzipEncoded=False,bundledSourceFiles=False))
 write(temporary/'manifest.json',dict(schemaVersion=1,records=records));os.rename(temporary,a.out)
 print(json.dumps(dict(status='retained-complete-expression-evidence',archiveBytes=archive.stat().st_size if archive.exists() else (a.out/'results.tar.gz').stat().st_size,members=len(members))))
if __name__=='__main__':main()
