#!/usr/bin/env python3
"""Archive complete receipts/checks; retain bulky native arrays externally by hash."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--ql-root',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args();sha=lambda b:hashlib.sha256(b).hexdigest()
if a.verify:
 m=json.loads((a.out/'manifest.json').read_text())
 for e in m['entries']:
  raw=(a.out/e['path']).read_bytes();assert sha(raw)==e['sha256']
  if 'logicalSHA256' in e:assert sha(gzip.decompress(raw))==e['logicalSHA256']
 assert {str(p.relative_to(a.out)) for p in a.out.rglob('*') if p.is_file()}=={e['path'] for e in m['entries']}|{'manifest.json'}
 print(json.dumps(dict(status='all-archive-members-verified',files=len(m['entries']))));raise SystemExit(0)
assert a.root and a.ql_root and not a.out.exists();a.out.mkdir(parents=True);entries=[];external=[]
def externalFile(path,role):
 raw=path.read_bytes();e=dict(path=str(path),bytes=len(raw),sha256=sha(raw),role=role)
 if path.suffix=='.gz' and not path.name.endswith('.tar.gz'):e['logicalSHA256']=sha(gzip.decompress(raw))
 external.append(e)
for path in sorted(a.root.rglob('*')):
 if not path.is_file():continue
 rel=path.relative_to(a.root)
 if rel.parts[0]=='edgeR':continue
 if path.name.endswith('.tar.gz'):
  externalFile(path,'Inspected reference package source cache; not a native dependency');continue
 if path.name=='native.json.gz' and rel.parts[0]!='work-sensitivity':
  externalFile(path,'Complete native moment/residual output');continue
 # First-pass per-gene checks are byte-identical to final outputs and retain
 # their own run receipts; avoid embedding the repeated detailed tables.
 if rel.parts[0] not in ['final','work-sensitivity'] and path.name in ['gene-checks.json.gz','independent-gene-errors.json.gz']:
  externalFile(path,'First-pass detailed checks');continue
 raw=path.read_bytes();compress=path.suffix in ['.log','.txt','.swift'] or (path.suffix=='.json' and len(raw)>131072)
 dest=a.out/(str(rel)+'.gz' if compress else str(rel));dest.parent.mkdir(parents=True,exist_ok=True);data=gzip.compress(raw,mtime=0) if compress else raw;dest.write_bytes(data)
 e=dict(path=str(dest.relative_to(a.out)),bytes=len(data),sha256=sha(data))
 if compress or path.suffix=='.gz':e['logicalSHA256']=sha(raw if compress else gzip.decompress(raw))
 entries.append(e)
protocol=json.loads((a.root/'final/protocol.json').read_text())
for c in protocol['cases']:
 for name in ['input.json.gz','reference.json.gz']:externalFile(a.ql_root/c['id']/name,'Previously qualified exact QL-stage input or fitted matrices')
for path,h in [(protocol['binary'],protocol['binarySHA256']),('/Users/n/numivivo-ql-moments-20260910/build/nb-moments',json.loads((a.root/'protocol.json').read_text())['binarySHA256']),('/Users/n/numivivo-ql-moments-20260910/sensitivity-build/nb-moments-work-sensitivity',json.loads((a.root/'work-sensitivity/protocol.json').read_text())['binarySHA256'])]:external.append(dict(path=path,sha256=h,role='Native compiled executable on '+protocol['host']))
(a.out/'manifest.json').write_text(json.dumps(dict(entries=entries,external=external,protocolSHA256=protocol['protocolSHA256'],qualification='Final per-gene checks, complete case receipts, independent checks, grid results, work sensitivity, underflow failure/repair and native logs archived. Complete bulk native arrays and prior QL fitted matrices remain externally retained at listed hashes.'),sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(files=len(entries),storedBytes=sum(e['bytes'] for e in entries),externalFiles=len(external))))
