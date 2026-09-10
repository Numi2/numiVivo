#!/usr/bin/env python3
"""Materialize a source-deduplicated archive as a regular native bundle atomically."""
import argparse,gzip,hashlib,json,os,shutil,tempfile
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();assert not a.out.exists(),'Output already exists'
m=json.loads((a.archive/'archive.json').read_text())
assert m['status']=='native-publication-and-replay-verified-before-archival'
assert m['source']=='../../annotated.h5ad'
limits={'plan.json':2097152,'report.json':536870912,'receipt.json':65536}
assert len(m['entries'])==3 and {e['logicalPath'] for e in m['entries']}==set(limits)
source=a.archive/m['source'];assert sha(source)==m['sourceSHA256'],'Shared source hash changed'
payloads={}
for e in m['entries']:
 assert e['path']==e['logicalPath']+'.gz'
 f=a.archive/e['path'];assert sha(f)==e['sha256'],'Compressed payload hash changed'
 with gzip.open(f,'rb') as stream:raw=stream.read(limits[e['logicalPath']]+1)
 assert len(raw)<=limits[e['logicalPath']] and hashlib.sha256(raw).hexdigest()==e['logicalSHA256'],'Logical payload changed or exceeds limit'
 payloads[e['logicalPath']]=raw
a.out.parent.mkdir(parents=True,exist_ok=True)
staging=Path(tempfile.mkdtemp(prefix='.numivivo-restore-',dir=a.out.parent))
try:
 shutil.copyfile(source,staging/'original.h5ad')
 assert sha(staging/'original.h5ad')==m['sourceSHA256'],'Copied source hash changed'
 for name,raw in payloads.items():(staging/name).write_bytes(raw)
 assert not a.out.exists(),'Output appeared during restoration'
 os.rename(staging,a.out)
finally:
 if staging.exists():shutil.rmtree(staging)
print(json.dumps(dict(status='materialized-for-native-replay',sourceSHA256=m['sourceSHA256'],output=str(a.out))))
