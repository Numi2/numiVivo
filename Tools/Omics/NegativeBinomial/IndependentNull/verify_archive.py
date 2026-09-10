#!/usr/bin/env python3
"""Verify every archived stored hash and original decoded byte identity."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);a=p.parse_args()
root=a.directory.resolve();manifest_bytes=(root/'manifest.json').read_bytes();manifest=json.loads(manifest_bytes)
sha=lambda b:hashlib.sha256(b).hexdigest()
assert len(manifest['files'])==manifest['fileCount']
assert len({f['path'] for f in manifest['files']})==manifest['fileCount']
for f in manifest['files']:
    path=(root/f['path']).resolve();assert root in path.parents
    stored=path.read_bytes();assert len(stored)==f['storedBytes'] and sha(stored)==f['storedSHA256']
    raw=gzip.decompress(stored) if f['path'].endswith('.gz') and not f['sourcePath'].endswith('.gz') else stored
    assert len(raw)==f['sourceBytes'] and sha(raw)==f['sourceSHA256']
assert sum(f['storedBytes'] for f in manifest['files'])==manifest['storedBytes']
print(json.dumps(dict(status='verified',files=manifest['fileCount'],storedBytes=manifest['storedBytes'],manifestSHA256=sha(manifest_bytes)),indent=2))
