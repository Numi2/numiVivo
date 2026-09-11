"""Verify compact source records and optionally every retained external artifact."""
import argparse, gzip, hashlib, json
from pathlib import Path, PurePosixPath

p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);p.add_argument('--external-root',type=Path);a=p.parse_args()
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def safe(name):
 path=PurePosixPath(name)
 assert not path.is_absolute() and '..' not in path.parts and str(path)==name
 return name
m=json.loads((a.evidence/'manifest.json').read_bytes());members=0
assert {x.name for x in a.evidence.iterdir()}=={'manifest.json',*(r['path'] for r in m['records'])}
for record in m['records']:
 path=a.evidence/safe(record['path'])
 assert path.stat().st_size==record['bytes'] and digest(path)==record['SHA256']
 raw=gzip.decompress(path.read_bytes())
 assert len(raw)==record['decodedBytes'] and hashlib.sha256(raw).hexdigest()==record['decodedSHA256']
 source=json.loads(raw)['sourceFiles'];assert len(source)==record['members']
 assert len({x['path'] for x in source})==len(source)
 for item in source:
  safe(item['path']);data=item['rawUTF8'].encode()
  assert len(data)==item['bytes'] and hashlib.sha256(data).hexdigest()==item['SHA256']
 members+=len(source)
assert len({x['path'] for x in m['externalArtifacts']})==len(m['externalArtifacts'])
if a.external_root:
 for item in m['externalArtifacts']:
  path=a.external_root/safe(item['path'])
  assert path.stat().st_size==item['bytes'] and digest(path)==item['SHA256'],item['path']
print(json.dumps(dict(status='passed',compactGroups=len(m['records']),sourceFiles=members,
 externalFiles=len(m['externalArtifacts']),externalArtifactsVerified=a.external_root is not None)))
