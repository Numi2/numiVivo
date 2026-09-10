#!/usr/bin/env python3
"""Archive every stage-audit output and verify compact evidence hashes."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
owner=Path(__file__).resolve().parent
def sha(b):return hashlib.sha256(b).hexdigest()
declaration=json.loads((a.root/'declaration.json').read_text())
for file,key in [('PROTOCOL.md','protocolSHA256'),('run.py','runnerSHA256'),('stages.R','referenceSHA256')]:
 assert sha((owner/file).read_bytes())==declaration[key]
trend=json.loads((a.root/'trend-declaration.json').read_text())
for file,digest in trend.items():assert sha((owner/file).read_bytes())==digest
assert json.loads((a.root/'complete.json').read_text())['status']=='completed-all-twenty-stage-audits'
assert json.loads((a.root/'trend-decomposition/complete.json').read_text())['status']=='completed-all-eighty-curves'
a.out.mkdir(parents=True,exist_ok=False);entries=[]
for file in sorted(a.root.rglob('*')):
 if not file.is_file():continue
 raw=file.read_bytes();relative=file.relative_to(a.root).as_posix()
 name=relative if file.suffix=='.gz' else relative+'.gz'
 data=raw if file.suffix=='.gz' else gzip.compress(raw,mtime=0)
 target=a.out/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
 entries.append(dict(path=name,bytes=len(data),sha256=sha(data),logicalPath=relative,logicalBytes=len(raw),logicalSHA256=sha(raw)))
manifest=dict(status='complete-diagnostic-evidence',entries=entries,sourceCommit='6e6f8b20338600d7965ce00d0b7b2e23c7c9c5d0',
 qualification='Post-score stage diagnosis; native inference unchanged; original count and model evidence remains in NullBenchmark.')
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for entry in entries:
 data=(a.out/entry['path']).read_bytes();assert sha(data)==entry['sha256'] and len(data)==entry['bytes']
 raw=data if entry['path']==entry['logicalPath'] else gzip.decompress(data)
 assert sha(raw)==entry['logicalSHA256'] and len(raw)==entry['logicalBytes']
print(json.dumps(dict(status='archived-and-verified',files=len(entries),bytes=sum(e['bytes'] for e in entries))))
