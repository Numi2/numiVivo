#!/usr/bin/env python3
"""Validate fixture identities before any native-versus-oracle comparison."""
import hashlib
import json
import pathlib
import sys


def verify(directory):
    root = pathlib.Path(directory).resolve(strict=True)
    manifest = json.loads((root / 'manifest.json').read_text())
    if (manifest.get('schema') != 'numivivo.org/external-conformance/v1'
            or manifest.get('oracle') != 'PySCF' or manifest.get('oracleVersion') != '2.8.0'):
        raise ValueError('Unsupported or unpinned oracle identity')
    expected = {'h2.json', 'lih.json', 'water.json'}
    entries = manifest.get('fixtures', [])
    if len(entries) != len(expected) or {entry['file'] for entry in entries} != expected:
        raise ValueError('The complete three-molecule fixture set is required')
    for entry in entries:
        path = root / entry['file']
        if path.resolve().parent != root or not path.is_file():
            raise ValueError('Fixture path escapes the reference directory')
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != entry['sha256']:
            raise ValueError(f"Fixture integrity failure: {entry['file']}")
    return manifest


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: verify_reference_manifest.py reference-directory')
    verify(sys.argv[1])
    print('Verified PySCF 2.8.0 fixture set and SHA-256 identities')
