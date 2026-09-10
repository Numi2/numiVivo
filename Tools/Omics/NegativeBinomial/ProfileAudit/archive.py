#!/usr/bin/env python3
"""Preserve complete profile evidence, including failed attempts and fallback cases."""
import argparse,gzip,hashlib,json,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
owner=Path(__file__).resolve().parent
def sha(b):return hashlib.sha256(b).hexdigest()
for declaration in [a.root/'declaration.json',a.root/'sensitivity/declaration.json']:
 for name,digest in json.loads(declaration.read_text()).items():assert sha((owner/name).read_bytes())==digest
assert json.loads((a.root/'reference-complete.json').read_text())['status']=='completed-all-twenty-reference-profile-audits'
assert json.loads((a.root/'sensitivity/complete.json').read_text())['status']=='completed-all-twenty-optimization-sensitivities'
summary=json.loads((a.root/'summary.json').read_text())
assert len(summary['perSplit'])==20 and summary['nativeValidPoints']==2305
a.out.mkdir(parents=True,exist_ok=False);entries=[]
for file in sorted(a.root.rglob('*')):
 if not file.is_file():continue
 raw=file.read_bytes();relative=file.relative_to(a.root).as_posix();name=relative if file.suffix=='.gz' else relative+'.gz'
 target=a.out/name;target.parent.mkdir(parents=True,exist_ok=True)
 if file.suffix=='.gz' and sys.platform=='darwin':
  subprocess.run(['/bin/cp','-c',str(file),str(target)],check=True)
  assert file.stat().st_ino!=target.stat().st_ino
 else:target.write_bytes(raw if file.suffix=='.gz' else gzip.compress(raw,mtime=0))
 stored=target.read_bytes()
 entries.append(dict(path=name,bytes=len(stored),sha256=sha(stored),logicalPath=relative,logicalBytes=len(raw),logicalSHA256=sha(raw)))
manifest=dict(status=summary['status'],entries=entries,sourceCommit='0c141130df7bb91d86ee42ff88202b5c899568c3',
 qualification='All attempts retained, including nine parametric sensitivity failures. Native numerical checks do not qualify calibration or power. Original H5AD/count evidence remains in NullBenchmark.')
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for e in entries:
 stored=(a.out/e['path']).read_bytes();assert sha(stored)==e['sha256'] and len(stored)==e['bytes']
 raw=stored if e['path']==e['logicalPath'] else gzip.decompress(stored)
 assert sha(raw)==e['logicalSHA256'] and len(raw)==e['logicalBytes']
print(json.dumps(dict(status='archived-and-verified',scientificStatus=summary['status'],files=len(entries),logicalStoredBytes=sum(e['bytes'] for e in entries))))
