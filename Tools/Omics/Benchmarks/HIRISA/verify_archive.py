#!/usr/bin/env python3
"""Check every archived member's stored and decoded bytes."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    root = args.archive
    manifest = json.loads((root / 'manifest.json').read_text())
    seen = set()
    for record in manifest['records']:
        name = record['storedPath']
        assert name not in seen
        seen.add(name)
        path = root / name
        assert path.resolve().is_relative_to(root.resolve()) and not path.is_symlink()
        stored = path.read_bytes()
        assert len(stored) == record['storedBytes']
        assert hashlib.sha256(stored).hexdigest() == record['storedSHA256']
        decoded = gzip.decompress(stored) if record['gzipEncoded'] else stored
        assert len(decoded) == record['sourceBytes']
        assert hashlib.sha256(decoded).hexdigest() == record['sourceSHA256']
    actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    assert actual == seen | {'manifest.json'}
    print(json.dumps({'status': 'passed', 'members': len(seen)}))

if __name__ == '__main__':
    main()
