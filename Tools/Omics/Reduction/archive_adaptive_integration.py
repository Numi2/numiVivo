#!/usr/bin/env python3
"""Retain exact evidence with deterministic compression and content deduplication.

Each source is LABEL=DIRECTORY. Original H5AD inputs and executables are external;
their canonical locations/hashes must be supplied in the external-inputs JSON.
Verification decompresses every stored object and checks both sets of hashes.
"""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--source',action='append',default=[])
p.add_argument('--external',type=Path)
p.add_argument('--verify',action='store_true')
p.add_argument('--reuse-gzip-from',type=Path,help='Reuse byte-identical raw evidence from an earlier archive, preserving its gzip bytes for Git deduplication.')
a=p.parse_args()
def sha(data):return hashlib.sha256(data).hexdigest()
if not a.verify:
 assert a.source and a.external
 a.out.mkdir(parents=True,exist_ok=False);(a.out/'objects').mkdir();entries=[];seen=set()
 for item in a.source:
  label,source=item.split('=',1);root=Path(source)
  assert label and '/' not in label and label not in ['.','..'] and label not in seen;seen.add(label)
  assert root.is_dir()
  for path in sorted(root.rglob('*')):
   if not path.is_file() or path.name.endswith('.h5ad') or '__pycache__' in path.parts:continue
   assert not path.is_symlink()
   raw=path.read_bytes();digest=sha(raw);object_path='objects/'+digest+'.gz';target=a.out/object_path
   if not target.exists():target.write_bytes(gzip.compress(raw,compresslevel=6,mtime=0))
   stored=target.read_bytes();entries.append(dict(path=label+'/'+path.relative_to(root).as_posix(),bytes=len(raw),sha256=digest,object=object_path,storedBytes=len(stored),storedSHA256=sha(stored)))
 manifest=dict(schemaVersion=1,encoding='gzip-content-addressed',files=entries,externalInputs=json.loads(a.external.read_text()),qualification='Exact evidence bytes with deduplicated storage. External source H5AD files and executable are not bundled. Native replay requires their restoration and the recorded executable identity.')
 (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
manifest=json.loads((a.out/'manifest.json').read_text());verified={};paths=set()
if a.reuse_gzip_from:
 prior_path=a.reuse_gzip_from/'archive-manifest.json';prior_bytes=prior_path.read_bytes()
 prior={e['uncompressedSHA256']:e for e in json.loads(prior_bytes) if e['encoding']=='gzip'};reused={}
 for entry in manifest['files']:
  digest=entry['sha256']
  if digest not in prior:continue
  if digest not in reused:
   old=prior[digest];relative=Path(old['path']);assert not relative.is_absolute() and '..' not in relative.parts
   stored=(a.reuse_gzip_from/relative).read_bytes();raw=gzip.decompress(stored)
   assert len(stored)==old['storedBytes'] and sha(stored)==old['storedSHA256']
   assert len(raw)==entry['bytes'] and sha(raw)==digest
   assert entry['object']=='objects/'+digest+'.gz'
   (a.out/entry['object']).write_bytes(stored);reused[digest]=(len(stored),sha(stored))
  entry['storedBytes'],entry['storedSHA256']=reused[digest]
 manifest['gzipReuseSource']=dict(archive=a.reuse_gzip_from.name,manifestSHA256=sha(prior_bytes),objects=len(reused),qualification='Each reused gzip stream independently decompresses to the exact recorded raw bytes. No external archive dependency is introduced.')
 (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
for entry in manifest['files']:
 assert entry['path'] not in paths;paths.add(entry['path'])
 assert entry['object']=='objects/'+entry['sha256']+'.gz'
 if entry['object'] not in verified:
  stored=(a.out/entry['object']).read_bytes();raw=gzip.decompress(stored)
  verified[entry['object']]=(len(raw),sha(raw),len(stored),sha(stored))
 assert verified[entry['object']]==(entry['bytes'],entry['sha256'],entry['storedBytes'],entry['storedSHA256'])
for entry in manifest.get('humanReadableFiles',[]):
 relative=Path(entry['path']);assert not relative.is_absolute() and '..' not in relative.parts
 raw=(a.out/relative).read_bytes();assert len(raw)==entry['bytes'] and sha(raw)==entry['sha256']
result=dict(status='passed',files=len(paths),uniqueObjects=len(verified),logicalBytes=sum(f['bytes'] for f in manifest['files']),storedBytes=sum(v[2] for v in verified.values()),humanReadableFiles=len(manifest.get('humanReadableFiles',[])),allRawAndStoredHashesVerified=True)
(a.out/'archive-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
