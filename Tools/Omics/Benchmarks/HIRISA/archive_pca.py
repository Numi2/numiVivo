#!/usr/bin/env python3
"""Archive full HIRISA PCA evidence with bounded, independently checked chunks."""
import argparse, gzip, hashlib, json
from pathlib import Path
from acquire import digest

def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();base=a.root/'pca';native=base/'quality-repair';out=a.output
 comparison=read(base/'comparison.json');assert comparison['status']=='passed' and comparison['cells']==1612594
 assert read(native/'full-publish-status.json')['returnCode']==read(native/'full-verify-status.json')['returnCode']==0
 assert '36 tests in 12 suites passed' in (base/'scale-tests-quality.log').read_text()
 for name in ('baron','hagai'):
  assert read(native/('regression-'+name)/'checks.json')['status']=='passed'
  assert read(native/('regression-'+name+'-exact.json'))['status']=='passed'
 assert read(native/'fixture-gate/checks.json')['status']=='passed'
 freeze=read(native/'full-output-freeze.json')
 for name,record in freeze.items():
  if name!='original.h5ad':assert digest(native/'full'/name)==record['SHA256']
 private_environments={}
 for source in sorted(base.rglob('*.json')):
  if source.name not in ('environment.json','quality-runtime-environment.json') or 'public-runtime-manifests' in source.parts:continue
  value=read(source);hardware=value.get('hardware')
  if not isinstance(hardware,str):continue
  labels=('Serial Number (system)','Hardware UUID','Provisioning UDID')
  redacted=[line.split(':',1)[0].strip() for line in hardware.splitlines() if ':' in line and line.split(':',1)[0].strip() in labels]
  if not redacted:continue
  value['hardware']='\n'.join(line.split(':',1)[0]+': [redacted]' if ':' in line and line.split(':',1)[0].strip() in labels else line for line in hardware.splitlines())+'\n'
  original=str(source.relative_to(base));value['originalPrivateManifestSHA256']=digest(source)
  target=base/'public-runtime-manifests'/original;target.parent.mkdir(parents=True,exist_ok=True)
  content=json.dumps(value,sort_keys=True,indent=2)+'\n'
  if target.exists():assert target.read_text()==content
  else:target.write_text(content)
  private_environments[original]=dict(SHA256=digest(source),publishedPath=str(target.relative_to(base)),redactedFields=redacted)
 paths=[]
 for p in sorted(base.rglob('*')):
  if not p.is_file():continue
  relative=p.relative_to(base)
  if str(relative) in private_environments:continue
  if any(part.startswith('.numivivo-') or part=='__pycache__' for part in relative.parts):continue
  if p.name=='numivivo' or p.suffix in ('.h5ad','.npy'):continue
  assert not p.is_symlink(),p
  paths.append(p)
 out.mkdir(parents=True,exist_ok=False);records=[];sources={};chunk_bytes=64*1024*1024
 for source in paths:
  name=str(source.relative_to(base));size=source.stat().st_size;before=digest(source);sources[name]=dict(bytes=size,SHA256=before)
  with source.open('rb') as src:
   offset=0;part=0
   while offset<size or (size==0 and part==0):
    length=min(chunk_bytes,size-offset);stored=name+('.part-%04d'%part if size>chunk_bytes else '')+'.gz'
    destination=out/stored;destination.parent.mkdir(parents=True,exist_ok=True);h=hashlib.sha256();remaining=length
    with destination.open('xb') as dst,gzip.GzipFile(fileobj=dst,mode='wb',filename='',mtime=0) as zipped:
     while remaining:
      block=src.read(min(1048576,remaining));assert block;zipped.write(block);h.update(block);remaining-=len(block)
    records.append(dict(sourcePath=name,sourceOffset=offset,sourceBytes=length,sourceSHA256=h.hexdigest(),
      storedPath=stored,storedBytes=destination.stat().st_size,storedSHA256=digest(destination),gzipEncoded=True))
    offset+=length;part+=1
   assert not src.read(1)
  assert source.stat().st_size==size and digest(source)==before
 manifest=dict(schemaVersion=2,records=records,fullSources=sources,sourceSHA256=comparison['sourceSHA256'],
  archiveScriptSHA256=digest(Path(__file__)),nativeReplayPassed=True,independentComparison=comparison,privateEnvironmentManifests=private_environments,
  externalPayloads='Original H5AD, native executables and reconstructible reference CSR .npy files remain external; their hashes are retained. Public runtime manifests omit device serial numbers and hardware/provisioning identifiers, with original private manifest hashes retained.',
  scope='Complete-source PCA numerical qualification and bounded storage admission; all failed and interrupted attempts retained. No million-cell graph/integration, biological validity or controlled CPU/GPU speed comparison.')
 (out/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
 print(json.dumps(dict(members=len(records),files=len(sources),storedBytes=sum(x['storedBytes'] for x in records),manifestSHA256=digest(out/'manifest.json'))))

if __name__=='__main__':main()
