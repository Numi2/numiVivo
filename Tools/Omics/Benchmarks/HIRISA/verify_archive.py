#!/usr/bin/env python3
"""Check every archived member's stored and decoded bytes."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def identity(stream):
    count = 0;sha = hashlib.sha256()
    for block in iter(lambda: stream.read(1_048_576), b''):
        count += len(block);sha.update(block)
    return count, sha.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    root = args.archive
    manifest = json.loads((root / 'manifest.json').read_text())
    seen = set()
    source_parts = {}
    for record in manifest['records']:
        name = record['storedPath']
        assert name not in seen
        seen.add(name)
        path = root / name
        assert path.resolve().is_relative_to(root.resolve()) and not path.is_symlink()
        with path.open('rb') as stream:
            assert identity(stream) == (record['storedBytes'], record['storedSHA256'])
        opener = gzip.open if record['gzipEncoded'] else open
        with opener(path, 'rb') as stream:
            assert identity(stream) == (record['sourceBytes'], record['sourceSHA256'])
        source_parts.setdefault(record['sourcePath'], []).append(record)
    # Version 2 can split large evidence files into bounded Git-friendly chunks.
    # Verify exact contiguous reconstruction as well as each stored/decoded part.
    if 'fullSources' in manifest:
        assert set(source_parts) == set(manifest['fullSources'])
        for name, expected in manifest['fullSources'].items():
            count = 0; sha = hashlib.sha256()
            for part in sorted(source_parts[name], key=lambda p: p['sourceOffset']):
                assert part['sourceOffset'] == count
                opener = gzip.open if part['gzipEncoded'] else open
                with opener(root / part['storedPath'], 'rb') as stream:
                    for block in iter(lambda: stream.read(1_048_576), b''):
                        count += len(block); sha.update(block)
            assert (count, sha.hexdigest()) == (expected['bytes'], expected['SHA256'])
    actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    assert actual == seen | {'manifest.json'}
    print(json.dumps({'status': 'passed', 'members': len(seen)}))

if __name__ == '__main__':
    main()
