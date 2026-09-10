#!/usr/bin/env python3
"""Archive reproducible pseudobulk inputs, all QL tables and complete check receipts."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args();sha=lambda b:hashlib.sha256(b).hexdigest()
if a.verify:
 m=json.loads((a.out/'manifest.json').read_text())
 for e in m['entries']:
  raw=(a.out/e['path']).read_bytes();assert sha(raw)==e['sha256']
  if 'logicalSHA256' in e:assert sha(gzip.decompress(raw))==e['logicalSHA256']
 assert {str(p.relative_to(a.out)) for p in a.out.rglob('*') if p.is_file()}=={e['path'] for e in m['entries']}|{'manifest.json'}
 print(json.dumps(dict(status='all-archive-members-verified',files=len(m['entries']))));raise SystemExit(0)
assert a.root and not a.out.exists();a.out.mkdir(parents=True);entries=[];external=[]
def write(rel,raw,compress=False):
 dest=a.out/(str(rel)+'.gz' if compress else str(rel));dest.parent.mkdir(parents=True,exist_ok=True);data=gzip.compress(raw,mtime=0) if compress else raw;dest.write_bytes(data)
 e=dict(path=str(dest.relative_to(a.out)),sha256=sha(data),bytes=len(data))
 if compress or dest.suffix=='.gz':e['logicalSHA256']=sha(raw if compress else gzip.decompress(data))
 entries.append(e)
for path in sorted(a.root.rglob('*')):
 if not path.is_file():continue
 rel=path.relative_to(a.root);raw=path.read_bytes()
 if path.name in ['reference.json.gz','native.json.gz']:
  external.append(dict(path=str(path),bytes=len(raw),sha256=sha(raw),logicalSHA256=sha(gzip.decompress(raw)),role='Complete reference matrices or native per-observation outputs'))
  if path.name=='reference.json.gz':
   result=json.loads(gzip.decompress(raw))
   for v in result['results'].values():
    for key in ['means','unitDF','unitDeviance','leverage','coefficients']:v.pop(key,None)
   write(rel.with_name('reference-tables.json'),json.dumps(result,sort_keys=True,separators=(',',':')).encode(),True)
  continue
 write(rel,raw,path.suffix in ['.log','.txt'] or (path.suffix=='.json' and len(raw)>131072))
protocol=json.loads((a.root/'protocol.json').read_text());native=json.loads((a.root/'native-summary.json').read_text())
external.extend(dict(path=c['source'],sha256=c['sourceSHA256'],logicalSHA256=c['logicalSHA256'],role='Original source counts/model report on '+protocol['host']) for c in protocol['cases'])
external.append(dict(path='/Users/n/numivivo-quasi-likelihood-20260910/build/nb-deviance',sha256=native['binarySHA256'],role='Native qualified primitive executable on '+protocol['host']))
(a.out/'manifest.json').write_text(json.dumps(dict(entries=entries,external=external,protocolSHA256=protocol['protocolSHA256'],qualification='All pseudobulk inputs, QL tables, source checks, arithmetic checks, native check summaries, boundary outputs and native logs archived; complete fitted matrices and per-observation native arrays retained externally at listed hashes'),sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(files=len(entries),storedBytes=sum(e['bytes'] for e in entries),externalFiles=len(external))))
