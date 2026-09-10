#!/usr/bin/env python3
"""Archive complete inputs/checks and bind external native model bytes."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if a.verify:
 m=json.loads((a.out/'manifest.json').read_text())
 for e in m['entries']:
  path=a.out/e['path'];assert sha(path)==e['sha256']
  if 'logicalSHA256' in e:assert hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest()==e['logicalSHA256']
 assert {str(p.relative_to(a.out)) for p in a.out.rglob('*') if p.is_file()}=={e['path'] for e in m['entries']}|{'manifest.json'}
 print(json.dumps(dict(status='all-archive-members-verified',files=len(m['entries']))));raise SystemExit(0)
assert a.root and not a.out.exists();a.out.mkdir(parents=True);entries=[];external=[]
protocol=json.loads((a.root/'protocol.json').read_text());runs=json.loads((a.root/'native-runs.json').read_text())
for rec in protocol['cases']:
 run=next(x for x in runs['runs'] if x['case']==rec['id'] and x['stage']=='fit')
 external.append(dict(path=str(Path(run['command'][-1]).with_name('model.json.gz')),sha256=run['outputSHA256'],logicalSHA256=run['logicalSHA256'],logicalBytes=run['outputBytes'],role='Complete frozen native model; verified through reference check'))
external.append(dict(path=runs['runs'][0]['command'][2],sha256=runs['binarySHA256'],role='Compiled native measurement harness'))
for path in sorted(a.root.rglob('*')):
 if not path.is_file():continue
 rel=path.relative_to(a.root)
 if rel.parts[0]=='build' or path.name in {'model.json.gz','probe-model.json.gz'}:
  external.append(dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size));continue
 raw=path.read_bytes();compress=path.suffix in {'.log','.txt','.tsv'} or path.suffix=='.json' and len(raw)>131072
 target=a.out/(str(rel)+'.gz' if compress else rel);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(gzip.compress(raw,mtime=0) if compress else raw)
 e=dict(path=str(target.relative_to(a.out)),bytes=target.stat().st_size,sha256=sha(target))
 if compress or path.suffix=='.gz':e['logicalSHA256']=hashlib.sha256(raw if compress else gzip.decompress(raw)).hexdigest()
 entries.append(e)
for name in ['sourceReport','productReport']:
 path=Path(protocol[name]);assert sha(path)==protocol[name+'SHA256'];external.append(dict(path=str(path),sha256=sha(path),role=name))
(a.out/'manifest.json').write_text(json.dumps(dict(entries=entries,external=external,binarySHA256=runs['binarySHA256'],protocolSHA256=protocol['protocolSHA256'],qualification='Exact training/test count inputs, native scores and independent checks; complete fitted models and original source reports externally retained.'),sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(files=len(entries),bytes=sum(e['bytes'] for e in entries),externalFiles=len(external))))
