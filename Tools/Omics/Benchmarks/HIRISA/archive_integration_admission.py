#!/usr/bin/env python3
"""Archive integration witness admission and lifecycle evidence after all checks."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1_048_576), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root, args.output
    read = lambda p: json.loads((root / p).read_text())
    native = read('admission-validation-complete.json')
    fixtures = read('admission-fixture-complete.json')
    assert native['status'] == fixtures['status'] == 'passed'
    assert native['binarySHA256'] == fixtures['binarySHA256']
    assert all(c['returnCode'] == 0 for c in native['commands'] + fixtures['commands'])
    assert '13 tests in 4 suites passed' in (root / 'admission-tests.log').read_text()
    assert "Build of product 'numivivo' complete!" in (root / 'admission-release-build.log').read_text()
    for name in ['admission-fixtures', 'mnn-fixtures']:
        result = read(name + '/checks.json')
        assert result['status'] == 'passed' and result['binarySHA256'] == native['binarySHA256']
    environment = read('runtime-release/environment.json')
    assert environment['binarySHA256'] == native['binarySHA256']
    assert not any(k in environment['hardware'] for k in ['Serial Number', 'Hardware UUID', 'Provisioning UDID'])
    selected = []
    for path in sorted(root.rglob('*')):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if '__pycache__' in relative.parts or path.name == 'numivivo':
            continue
        assert not path.is_symlink() and path.stat().st_size <= 64 * 1024 * 1024, path
        selected.append(path)
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for source in selected:
        name = str(source.relative_to(root))
        raw = source.read_bytes()
        stored = name + '.gz'
        target = output / stored
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as f:
            with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0, compresslevel=6) as z:
                z.write(raw)
        source_sha = hashlib.sha256(raw).hexdigest()
        assert digest(source) == source_sha
        records.append(dict(sourcePath=name, sourceBytes=len(raw), sourceSHA256=source_sha,
                            storedPath=stored, storedBytes=target.stat().st_size,
                            storedSHA256=digest(target), gzipEncoded=True))
    manifest = dict(schemaVersion=1, records=records, archiveScriptSHA256=digest(Path(__file__)),
                    binarySHA256=native['binarySHA256'],
                    scope='Shared integration witness admission, native solver regression tests, release build, fitted/query ridge and MNN lifecycle, independent numerical reference and corruption controls. Historical full-HIRISA admission audit is retained; no full-HIRISA integration execution or biological qualification is claimed.')
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(json.dumps(dict(members=len(records), storedBytes=sum(r['storedBytes'] for r in records),
                         manifestSHA256=digest(output / 'manifest.json'))))


if __name__ == '__main__':
    main()
