#!/usr/bin/env python3
"""Archive exact GO evidence, splitting large files into reversible bounded chunks."""
import argparse,gzip,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work',type=Path);p.add_argument('--out',type=Path,required=True)
p.add_argument('--verify',action='store_true')
a=p.parse_args();base=Path(__file__).parents[2]/'Reduction/archive_adaptive_integration.py'
CHUNK=64*1024*1024
if not a.verify:
 if a.work is None:p.error('packing requires --work')
 with tempfile.TemporaryDirectory(prefix='numivivo-go-archive-') as directory:
  stage=Path(directory);root=stage/'run';root.mkdir();large=[]
  selected=['annotations','source','prepared','prepared-repeat','first','repeat','scores','scores-repeat',
            'fixed-control','fixed-control-repeat','fixed-control-scores','fixed-control-scores-repeat']
  files=[f for label in selected for f in sorted((a.work/label).rglob('*')) if f.is_file()]
  files += [f for f in sorted(a.work.iterdir()) if f.is_file() and f.suffix in ('.json','.md','.log') and not f.name.startswith('packaging')]
  for source in files:
   assert not source.is_symlink();relative=source.relative_to(a.work)
   if source.stat().st_size<=CHUNK:
    dest=root/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest);continue
   dest=root/relative.parent/(relative.name+'.parts');dest.mkdir(parents=True)
   whole=hashlib.sha256();parts=[];size=0
   with source.open('rb') as stream:
    while raw:=stream.read(CHUNK):
     path=dest/(f'{len(parts):04d}.part');path.write_bytes(raw);whole.update(raw);size+=len(raw)
     parts.append('run/'+path.relative_to(root).as_posix())
   large.append(dict(path='run/'+relative.as_posix(),bytes=size,sha256=whole.hexdigest(),parts=parts))
  external=stage/'external.json';external.write_text(json.dumps(dict(referencePath='/Users/home/numivivo-norman-prepared-20260909/reference.npz',referenceSHA256='ec22f61196cf51f7c1a62c681e4076bd9ea269429f1423dfd65e21c11a7ce8bf',baseCommit='61dcea919e3cf2862bfc4910e37b060520fde928')))
  subprocess.run([sys.executable,str(base),'--out',str(a.out),'--external',str(external),'--source','run='+str(root)],check=True)
  manifest=json.loads((a.out/'manifest.json').read_text());manifest['largeArtifacts']=large
  manifest['largeArtifactEncoding']='Concatenate decompressed parts in order; verify the recorded whole-file hash and size. Original bytes are not reserialized.'
  (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
subprocess.run([sys.executable,str(base),'--verify','--out',str(a.out)],check=True)
m=json.loads((a.out/'manifest.json').read_text());files={e['path']:e for e in m['files']}
for entry in m.get('largeArtifacts',[]):
 whole=hashlib.sha256();size=0
 for part in entry['parts']:
  e=files[part]
  with gzip.open(a.out/e['object'],'rb') as stream:
   while raw:=stream.read(1048576):whole.update(raw);size+=len(raw)
 assert size==entry['bytes'] and whole.hexdigest()==entry['sha256'],entry['path']
result=dict(status='passed',largeFilesReassembledAndHashVerified=len(m.get('largeArtifacts',[])),maximumStoredObjectBytes=max(e['storedBytes'] for e in m['files']))
assert result['maximumStoredObjectBytes']<100*1024*1024
(a.out/'chunk-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
