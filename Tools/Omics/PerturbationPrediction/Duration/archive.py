#!/usr/bin/env python3
"""Archive exact native snapshots and reports once by content; verify before restore."""
import argparse,gzip,hashlib,io,json,shutil,subprocess,tarfile
from pathlib import Path,PurePosixPath
from common import sha,write

def verify(archive):
 with tarfile.open(archive,'r:gz') as t:
  items=t.getmembers();assert len({i.name for i in items})==len(items) and all(i.isfile() for i in items)
  members=json.load(t.extractfile('members.json'));objects={v['SHA256'] for v in members.values()}
  assert {i.name for i in items}=={'members.json'}|{'objects/'+h for h in objects}
  for name,v in members.items():
   p=PurePosixPath(name);h=v['SHA256'];assert 0<=v['mode']<=0o777;assert not p.is_absolute() and '..' not in p.parts and str(p)==name
   assert len(h)==64 and all(c in '0123456789abcdef' for c in h) and t.getmember('objects/'+h).size==v['bytes']
  for item in items:
   if not item.name.startswith('objects/'):continue
   digest=hashlib.sha256()
   with t.extractfile(item) as f:
    for b in iter(lambda:f.read(1048576),b''):digest.update(b)
   assert digest.hexdigest()==item.name.split('/')[1]
 return members

def restore(archive,out):
 members=verify(archive);out.mkdir(parents=True,exist_ok=False);byhash={}
 for name,v in members.items():byhash.setdefault(v['SHA256'],[]).append(name)
 with tarfile.open(archive,'r:gz') as t:
  for item in t.getmembers():
   if not item.name.startswith('objects/'):continue
   h=item.name.split('/')[1];previous=None
   for name in byhash[h]:
    p=out/name;p.parent.mkdir(parents=True,exist_ok=True);cloned=False
    if previous is not None:cloned=subprocess.run(['/bin/cp','-c',str(previous),str(p)],capture_output=True).returncode==0
    if not cloned:
     assert not p.exists() and shutil.disk_usage(out).free>item.size+200000000
     with t.extractfile(item) as f,p.open('xb') as g:shutil.copyfileobj(f,g,1048576)
    assert sha(p)==h;p.chmod(members[name]['mode']);previous=p
 return members

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=['pack','verify','restore']);p.add_argument('--study',type=Path);p.add_argument('--out',type=Path);p.add_argument('--archive',type=Path);a=p.parse_args()
 if a.action=='verify':print(json.dumps(dict(status='verified',members=len(verify(a.archive)))));return
 if a.action=='restore':print(json.dumps(dict(status='restored',members=len(restore(a.archive,a.out)))));return
 a.out.mkdir(parents=True,exist_ok=False);files=[]
 for folder in ['inputs','native','scores','scores-repeat','failure-checks','plots']:
  files += [f for f in (a.study/folder).rglob('*') if f.is_file()]
 files += [f for f in a.study.iterdir() if f.is_file() and f.suffix in {'.json','.log'}]
 files += [a.study/'build/numivivo-omics',a.study/'build/sources.sha256']
 files=sorted(set(files));assert all(not f.is_symlink() for f in files)
 members={str(f.relative_to(a.study)):dict(SHA256=sha(f),bytes=f.stat().st_size,mode=f.stat().st_mode & 0o777) for f in files};objects={}
 for f in files:objects.setdefault(members[str(f.relative_to(a.study))]['SHA256'],f)
 archive=a.out/'results.tar.gz';raw=(json.dumps(members,sort_keys=True,indent=2)+'\n').encode()
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  info=tarfile.TarInfo('members.json');info.size=len(raw);t.addfile(info,io.BytesIO(raw))
  for h,f in sorted(objects.items()):
   info=tarfile.TarInfo('objects/'+h);info.size=f.stat().st_size
   with f.open('rb') as source:t.addfile(info,source)
 assert verify(archive)==members
 write(a.out/'contents.json',dict(status='verified',members=members,logicalBytes=sum(v['bytes'] for v in members.values()),uniqueBytes=sum(f.stat().st_size for f in objects.values()),archiveSHA256=sha(archive),scope='Exact native H5AD snapshots, plans, reports, receipts, full aggregate inputs, all predictions/outcomes and executable; source-raw provenance binds the prior external experiment.'))
 # Preserve authored recipe and new native owners as exact text alongside the commit.
 repo=Path(__file__).resolve().parents[4];recipefiles=[f for f in Path(__file__).parent.iterdir() if f.is_file()]
 recipefiles += [repo/n for n in ['Sources/NumiVivoKit/Omics/VivoDurationPerturbation.swift','Sources/NumiVivoKit/Omics/VivoDurationPerturbationIO.swift','Sources/NumiVivoKit/Omics/VivoPerturbationIO.swift','Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift','Tests/NumiVivoIntegrationTests/SingleCellDurationTests.swift','Tests/NumiVivoIntegrationTests/SingleCellPerturbationTests.swift','Tools/Omics/H5AD/build.sh']]
 sources=[]
 for f in sorted(set(recipefiles)):
  raw=f.read_bytes();sources.append(dict(path=str(f.relative_to(repo)),bytes=len(raw),SHA256=sha(f),rawUTF8=raw.decode()))
 recipe=(json.dumps(dict(schemaVersion=1,sourceFiles=sources),sort_keys=True,separators=(',',':'))+'\n').encode();(a.out/'recipe.json.gz').write_bytes(gzip.compress(recipe,compresslevel=9,mtime=0))
 records=[]
 for f in sorted(a.out.iterdir()):
  data=f.read_bytes();compressed=f.name=='recipe.json.gz';decoded=gzip.decompress(data) if compressed else data
  records.append(dict(sourcePath=f.name.removesuffix('.gz') if compressed else f.name,sourceBytes=len(decoded),sourceSHA256=hashlib.sha256(decoded).hexdigest(),storedPath=f.name,storedBytes=len(data),storedSHA256=sha(f),gzipEncoded=compressed,bundledSourceFiles=compressed))
 write(a.out/'manifest.json',dict(schemaVersion=1,records=records));print(json.dumps(dict(members=len(members),logicalBytes=sum(v['bytes'] for v in members.values()),archiveBytes=archive.stat().st_size,archiveSHA256=sha(archive))))
if __name__=='__main__':main()
