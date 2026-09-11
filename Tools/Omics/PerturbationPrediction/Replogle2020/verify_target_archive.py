#!/usr/bin/env python3
"""Verify stored, decoded and bundled original evidence bytes."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    root = args.archive.resolve()
    manifest = json.loads((root/'manifest.json').read_text())
    seen, sources, bundled = set(), set(), 0
    for item in manifest['records']:
        name = item['storedPath']
        assert name not in seen
        seen.add(name)
        path = root/name
        assert path.resolve().is_relative_to(root) and not path.is_symlink()
        raw = path.read_bytes()
        assert len(raw) == item['storedBytes'] and hashlib.sha256(raw).hexdigest() == item['storedSHA256']
        decoded = gzip.decompress(raw) if item['gzipEncoded'] else raw
        assert len(decoded) == item['sourceBytes'] and hashlib.sha256(decoded).hexdigest() == item['sourceSHA256']
        if item.get('bundledSourceFiles'):
            document = json.loads(decoded)
            assert document['schemaVersion'] == 1
            for source in document['sourceFiles']:
                assert source['path'] not in sources
                sources.add(source['path'])
                original = source['rawUTF8'].encode('utf-8')
                assert len(original) == source['bytes'] and hashlib.sha256(original).hexdigest() == source['SHA256']
                bundled += 1
    assert {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()} == seen|{'manifest.json'}
    print(json.dumps(dict(status='passed', members=len(seen), bundledOriginalFiles=bundled)))


if __name__ == '__main__':
    main()
