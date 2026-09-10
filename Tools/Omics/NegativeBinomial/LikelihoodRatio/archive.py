#!/usr/bin/env python3
"""Archive complete compact measurements, checks and preserved failed attempts."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args()
sha=lambda b:hashlib.sha256(b).hexdigest()
if a.verify:
 m=json.loads((a.out/'manifest.json').read_text())
 for e in m['entries']:
  raw=(a.out/e['path']).read_bytes();assert sha(raw)==e['sha256']
  if 'logicalSHA256' in e:assert sha(gzip.decompress(raw))==e['logicalSHA256']
 assert {str(p.relative_to(a.out)) for p in a.out.rglob('*') if p.is_file()}=={e['path'] for e in m['entries']}|{'manifest.json'}
 print(json.dumps(dict(status='all-archive-members-verified',files=len(m['entries']))));raise SystemExit(0)
assert a.root and not a.out.exists();a.out.mkdir(parents=True);entries=[];external=[]
for path in sorted(a.root.rglob('*')):
 if not path.is_file():continue
 rel=path.relative_to(a.root)
 # The separately mirrored product is identical to product/ above it.
 if str(rel).startswith('remote-execution/product/'):continue
 if path.name=='native.json.gz' or (rel.parts[0]=='initial-reference-checks' and path.name=='gene-checks.json.gz'):
  raw=path.read_bytes();external.append(dict(path=str(path),bytes=len(raw),sha256=sha(raw),logicalSHA256=sha(gzip.decompress(raw)),role='Complete native measurement or original per-gene reference pass, retained externally'));continue
 raw=path.read_bytes();compress=path.suffix in {'.log','.txt','.tsv','.swift'} or (path.suffix=='.json' and len(raw)>131072)
 dest=a.out/(str(rel)+'.gz' if compress else str(rel));dest.parent.mkdir(parents=True,exist_ok=True)
 dest.write_bytes(gzip.compress(raw,mtime=0) if compress else raw)
 e=dict(path=str(dest.relative_to(a.out)),bytes=dest.stat().st_size,sha256=sha(dest.read_bytes()))
 if compress or path.suffix=='.gz':e['logicalSHA256']=sha(raw if compress else gzip.decompress(raw))
 entries.append(e)
protocol=json.loads((a.root/'protocol.json').read_text());native=json.loads((a.root/'native-complete.json').read_text());product=json.loads((a.root/'product/complete.json').read_text())
external.extend(dict(path=c['source'],sha256=c['sourceSHA256'],logicalSHA256=c['logicalSHA256'],role='Original qualified source model with complete counts and metadata') for c in protocol['cases'])
external.extend([dict(path='/Users/n/numivivo-likelihood-ratio-20260910/build/nb-lrt',sha256=native['binarySHA256'],role='Native compiled measurement harness'),dict(path='/Users/n/numivivo-likelihood-ratio-20260910/numivivo',sha256=product['binarySHA256'],role='Qualified full product executable')])
(a.out/'manifest.json').write_text(json.dumps(dict(entries=entries,external=external,protocolSHA256=protocol['protocolSHA256'],qualification='Complete final per-gene checks, all initial check summaries, reference tables, product reports/receipts and retained failed attempts; complete native outputs, initial per-gene tables, source model archives and executables remain externally retained'),sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(files=len(entries),storedBytes=sum(e['bytes'] for e in entries),externalFiles=len(external))))
