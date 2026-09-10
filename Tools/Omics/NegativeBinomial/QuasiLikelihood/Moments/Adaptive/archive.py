#!/usr/bin/env python3
"""Archive adaptive qualification receipts and compact results with bulk arrays external."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args();sha=lambda b:hashlib.sha256(b).hexdigest()
if a.verify:
 m=json.loads((a.out/'manifest.json').read_text())
 for e in m['entries']:
  b=(a.out/e['path']).read_bytes();assert sha(b)==e['sha256']
  if 'logicalSHA256' in e:assert sha(gzip.decompress(b))==e['logicalSHA256']
 assert {str(p.relative_to(a.out)) for p in a.out.rglob('*') if p.is_file()}=={e['path'] for e in m['entries']}|{'manifest.json'}
 print(json.dumps(dict(status='all-archive-members-verified',files=len(m['entries']))));raise SystemExit(0)
assert a.root and not a.out.exists();a.out.mkdir(parents=True);entries=[];external=[]
for path in sorted(a.root.rglob('*')):
 if not path.is_file():continue
 rel=path.relative_to(a.root);raw=path.read_bytes()
 if path.name=='native.json.gz':
  external.append(dict(path=str(path),bytes=len(raw),sha256=sha(raw),logicalSHA256=sha(gzip.decompress(raw)),role='Full adaptive family moment/residual output'));continue
 compress=path.suffix in ['.log','.txt','.swift'] or (path.suffix=='.json' and len(raw)>32768)
 data=gzip.compress(raw,mtime=0) if compress else raw;dest=a.out/(str(rel)+'.gz' if compress else str(rel));dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
 e=dict(path=str(dest.relative_to(a.out)),bytes=len(data),sha256=sha(data))
 if compress or path.suffix=='.gz':e['logicalSHA256']=sha(raw if compress else gzip.decompress(raw))
 entries.append(e)
execution=json.loads((a.root/'execution.json').read_text());external.extend(execution['externalExecutables'])
(a.out/'manifest.json').write_text(json.dumps(dict(entries=entries,external=external,protocolSHA256=execution['protocolSHA256'],qualification='Complete native adaptive family, prior direct comparisons, 124 independent recovered moments, extended and legacy grids, derivative reference, native logs and retained compilation/transport observations. Bulk arrays remain external; parent direct evidence retains original resource failures.'),sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(files=len(entries),storedBytes=sum(e['bytes'] for e in entries),externalFiles=len(external))))
