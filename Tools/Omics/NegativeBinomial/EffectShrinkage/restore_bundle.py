#!/usr/bin/env python3
"""Restore a verified treatment archive using an explicitly supplied source H5AD."""
import argparse,gzip,hashlib,json,os,shutil,tempfile
from pathlib import Path

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
p=argparse.ArgumentParser(description=__doc__)
for name in ['archive','source','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();assert not a.out.exists()
m=json.loads((a.archive/'archive.json').read_text())
assert m['status']=='published-and-replay-verified-before-archival'
assert sha(a.source)==m['sourceSHA256']
limits={'plan.json':2097152,'report.json':536870912,'receipt.json':65536}
assert len(m['entries'])==3 and {e['logicalPath'] for e in m['entries']}==set(limits)
payloads={}
for e in m['entries']:
 assert e['path']==e['logicalPath']+'.gz'
 path=a.archive/e['path'];assert sha(path)==e['sha256']
 with gzip.open(path,'rb') as f:raw=f.read(limits[e['logicalPath']]+1)
 assert len(raw)<=limits[e['logicalPath']] and hashlib.sha256(raw).hexdigest()==e['logicalSHA256']
 payloads[e['logicalPath']]=raw
a.out.parent.mkdir(exist_ok=True,parents=True);staging=Path(tempfile.mkdtemp(prefix='.numivivo-restore-',dir=a.out.parent))
try:
 shutil.copyfile(a.source,staging/'original.h5ad');assert sha(staging/'original.h5ad')==m['sourceSHA256']
 for name,raw in payloads.items():(staging/name).write_bytes(raw)
 assert not a.out.exists();os.rename(staging,a.out)
finally:
 if staging.exists():shutil.rmtree(staging)
print(json.dumps(dict(status='materialized-for-native-replay',sourceSHA256=m['sourceSHA256'],output=str(a.out))))
