#!/usr/bin/env python3
"""Archive the actual cohort run, retaining failure evidence and external hashes."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(b):return hashlib.sha256(b).hexdigest()
entries=[]
files=[a.work/x for x in ['plan.json','cohort.json','audit.json','reference.npz','check.json','verify.json',
    'native-receipt.json','native-time.log','verify-time.log','release-build.log','release-build-final.log',
    'duplicate-cleanup.log','attempts.json','author-cell-population.py','author-LICENSE.txt','streaming-regression.log']]
for directory in ['native','selection-checks','selection-checks-final','selection-checks-verified','streaming-regression']:
    files+=sorted(x for x in (a.work/directory).rglob('*') if x.is_file())
for source in files:
    assert source.stat().st_size<150_000_000,source
    raw=source.read_bytes();rel=source.relative_to(a.work).as_posix()
    stored=gzip.compress(raw,mtime=0) if source.suffix not in ['.npz','.gz'] else raw
    destination=rel+'.gz' if stored is not raw else rel
    dest=a.out/destination;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(stored)
    entries.append(dict(path=destination,bytes=len(stored),sha256=sha(stored),logicalPath=rel,logicalBytes=len(raw),logicalSHA256=sha(raw)))
external=[]
for source in [a.work/'numivivo',a.work.parent/'geo-restoration/annotated.h5ad']:
    with source.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
    external.append(dict(path=str(source),bytes=source.stat().st_size,sha256=digest))
manifest=dict(schemaVersion=1,status='verified-native-author-cohort-counts-prediction-pending',entries=entries,external=external,
    retainedRemoteSource='/Users/n/numivivo-adamson-target-kernel-20260910/cohort/native/original.h5ad')
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
# Verify stored bytes and their decoded payloads independently of write success.
for e in entries:
    b=(a.out/e['path']).read_bytes();assert len(b)==e['bytes'] and sha(b)==e['sha256']
    raw=gzip.decompress(b) if e['path']!=e['logicalPath'] else b
    assert len(raw)==e['logicalBytes'] and sha(raw)==e['logicalSHA256']
print(json.dumps(dict(status='archived-and-verified',files=len(entries),storedBytes=sum(e['bytes'] for e in entries))))
