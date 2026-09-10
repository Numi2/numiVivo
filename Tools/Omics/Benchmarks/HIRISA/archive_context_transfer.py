#!/usr/bin/env python3
"""Archive all preparation-transfer folds and checks with bounded, verified chunks."""
import argparse
import gzip
import hashlib
from pathlib import Path
from prepare_context_predictions import read, sha, write

LINEAGES = ('B', 'Mono', 'NK', 'CD4-T', 'CD8-T', 'other-T')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args(); root = a.root
    complete = read(root/'prediction-complete.json')
    freeze = read(root/'prediction-output-freeze.json')
    scores = read(root/'scores.json')
    assert complete['allNativeReplaysPassed'] and complete['folds'] == 120
    assert scores['foldCount'] == 120 and len(scores['comparisons']) == 12
    assert scores['outputFreezeSHA256'] == complete['outputFreezeSHA256'] == sha(root/'prediction-output-freeze.json')
    assert read(root/'native-complete.json')['allSixLineagesAndReplaysPassed']
    assert read(root/'scorer-regression.json')['folds'] == 79
    for name, digest in freeze['files'].items():
        assert sha(root/name) == digest, name
    paths = {p for p in root.iterdir() if p.is_file() and p.suffix in {'.json', '.log', '.py', '.bin'}}
    restoration = {}
    source_sha = read(root/'freeze.json')['sourceSHA256']
    for lineage in LINEAGES:
        native, prediction = root/(lineage+'-native'), root/(lineage+'-prediction')
        paths.update(p for p in native.iterdir() if p.suffix == '.json')
        paths.update(p for p in prediction.rglob('*.json') if 'source' not in p.relative_to(prediction).parts)
        # Snapshot payloads are APFS clones of the retained original source. Source
        # report/plan/receipt copies are represented once, with exact restoration.
        for snapshot in (native/'original.h5ad', prediction/'source/original.h5ad'):
            assert sha(snapshot) == source_sha
            restoration[str(snapshot.relative_to(root))] = dict(externalSource='original full HIRISA H5AD', SHA256=source_sha, bytes=snapshot.stat().st_size)
        for name in ('plan.json', 'report.json', 'receipt.json'):
            source = prediction/'source'/name
            assert sha(source) == sha(native/name)
            restoration[str(source.relative_to(root))] = dict(sourcePath='study/'+str((native/name).relative_to(root)), SHA256=sha(source))
    for folder in ('prediction-inputs', 'reference'):
        paths.update(p for p in (root/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    sources = [('study/'+str(p.relative_to(root)), p) for p in paths]
    sources.append(('tools/archive_context_transfer.py', Path(__file__)))
    a.out.mkdir(parents=True, exist_ok=False)
    records, full = [], {}
    for name, source in sorted(sources):
        assert not source.is_symlink()
        digest = sha(source); size = source.stat().st_size
        full[name] = dict(bytes=size, SHA256=digest)
        offset = 0
        with source.open('rb') as stream:
            part = 0
            while True:
                raw = stream.read(16*2**20)
                if not raw and (part or size): break
                stored = name+f'.part-{part:04d}.gz'
                encoded = gzip.compress(raw, compresslevel=6, mtime=0)
                target = a.out/stored; target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(encoded)
                records.append(dict(sourcePath=name, sourceOffset=offset, sourceBytes=len(raw),
                    sourceSHA256=hashlib.sha256(raw).hexdigest(), storedPath=stored, storedBytes=len(encoded),
                    storedSHA256=hashlib.sha256(encoded).hexdigest(), gzipEncoded=True))
                offset += len(raw); part += 1
                if not raw: break
        assert offset == size and sha(source) == digest
    write(a.out/'manifest.json', dict(schemaVersion=2, records=records, fullSources=full,
        restoration=restoration, previousPreparationArchive='../2026-09-10-context-transfer-preparation/manifest.json',
        folds=120, completedFolds=scores['completedFolds'], nativeReplaysPassed=True,
        sourceSHA256=source_sha, outputFreezeSHA256=sha(root/'prediction-output-freeze.json'),
        scoresSHA256=sha(root/'scores.json'), scope=scores['scope']))
    print(dict(members=len(records), storedBytes=sum(r['storedBytes'] for r in records), manifestSHA256=sha(a.out/'manifest.json')))


if __name__ == '__main__':
    main()
