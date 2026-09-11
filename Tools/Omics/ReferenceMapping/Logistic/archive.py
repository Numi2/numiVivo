#!/usr/bin/env python3
"""Archive compact qualification evidence; bind full native bundles by hash."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    groups, records = {}, []

    def add(group, path, label):
        raw = path.read_bytes()
        groups.setdefault(group, []).append(dict(
            path=label, bytes=len(raw), SHA256=digest(raw), rawUTF8=raw.decode()))

    def store(name, raw, bundled=False):
        encoded = gzip.compress(raw, compresslevel=9, mtime=0)
        stored = name + '.gz'
        (args.out / stored).write_bytes(encoded)
        records.append(dict(sourcePath=name, sourceBytes=len(raw),
            sourceSHA256=digest(raw), storedPath=stored,
            storedBytes=len(encoded), storedSHA256=digest(encoded),
            gzipEncoded=True, bundledSourceFiles=bundled))

    text_suffixes = {'.json', '.log', '.txt', '.sh', '.sha256', '.py', '.swift', '.md'}
    for path in sorted(args.study.iterdir()):
        if path.is_file() and path.suffix in text_suffixes:
            add('execution', path, 'study/' + path.name)
    for directory in ['inputs', 'inputs-qualified', 'native', 'native-qualified',
                      'regression-metadata', 'regression-qualified-metadata',
                      'regression-fixtures', 'scores-qualified', 'scores-repeat',
                      'recipe', 'recipe-qualified']:
        for path in sorted((args.study / directory).rglob('*')):
            if path.is_file() and not path.is_symlink() and path.suffix in text_suffixes:
                add('execution', path, 'study/' + str(path.relative_to(args.study)))
    for path in sorted(Path(__file__).resolve().parent.iterdir()):
        if path.is_file():
            add('recipe', path, 'repo/' + str(path.relative_to(args.repo.resolve())))
    for name in [
        'Sources/NumiVivoKit/Omics/VivoReferenceLogistic.swift',
        'Sources/NumiVivoKit/Omics/VivoSingleCellReference.swift',
        'Sources/NumiVivoKit/Omics/VivoSingleCellReferenceIO.swift',
        'Tools/Omics/H5AD/build.sh',
        'Tests/NumiVivoIntegrationTests/ReferenceLogisticTests.swift',
        'Tests/NumiVivoIntegrationTests/SingleCellReferenceTests.swift',
    ]:
        add('native-source', args.repo / name, 'repo/' + name)
    for name, files in sorted(groups.items()):
        raw = (json.dumps(dict(schemaVersion=1, sourceFiles=files),
                          sort_keys=True, separators=(',', ':')) + '\n').encode()
        store(name + '.json', raw, bundled=True)
    # Keep all per-cell probability comparisons, not just aggregate scores.
    for path in sorted((args.study / 'scores-qualified').glob('*.npz')):
        store(path.name, path.read_bytes())
    (args.out / 'manifest.json').write_text(json.dumps(
        dict(schemaVersion=1, records=records), indent=2) + '\n')
    print(json.dumps(dict(groups=len(groups), files=sum(map(len, groups.values())),
                          storedBytes=sum(item['storedBytes'] for item in records))))


if __name__ == '__main__':
    main()
