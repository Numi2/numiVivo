#!/usr/bin/env python3
"""Archive qualified graph-reference reuse without duplicating parent payloads."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    study = args.root
    root = study/'optimized-clustering'
    qualification = read(root/'checker-qualification-v2/qualification.json')
    collection = read(root/'collector-baron/collection.json')
    assert qualification['status'] == 'passed'
    assert len(qualification['commands']) == 9
    assert sum(not v['expectedSuccess'] for v in qualification['commands']) == 7
    for name, digest in qualification['sourceIdentities'].items():
        assert sha(Path(name)) == digest, name
    assert read(root/'collector-qualification.json')['status'] == 'passed'
    assert read(root/'collector-baron/complete.json')['status'] == 'passed'
    assert read(root/'collector-baron/remote-before.json') == read(root/'collector-baron/remote-after.json') == collection['remotePayloads']
    for name, expected in collection['remotePayloads'].items():
        path = root/'collector-baron/bundle'/name
        assert dict(bytes=path.stat().st_size, SHA256=sha(path)) == expected
    assert read(root/'checker-qualification/attempt-status.json')['status'] == 'failed-test-harness-expectation'
    assert not (root/'collector-cross-cohort-rejected/bundle').exists()
    full_check = read(root/'full-reference-check/checks.json')
    full_reference = read(study/'clustering/independent-reference/reference.json')
    assert full_check['status'] == 'passed' and full_check['cells'] == 1612594
    assert full_check['referenceSHA256'] == sha(study/'clustering/independent-reference/reference.json')
    assert full_check['bindings']['partitionsSHA256'] == sha(study/'clustering/independent-reference/partitions.npz')
    assert full_reference['graphReceiptSHA256'] == sha(study/'graph/full/receipt.json')
    selected = []
    for folder in ['checker-qualification', 'checker-qualification-v2', 'collector-baron',
                   'collector-cross-cohort-rejected']:
        selected += [path for path in (root/folder).rglob('*') if path.is_file()
                     and 'bundle' not in path.relative_to(root/folder).parts]
    selected += [root/'collector-baron/bundle'/name for name in collection['copied']]
    selected += [root/name for name in ['collector-qualification.json',
                                       'collector-cross-cohort-rejected.log', 'continue_independent.py']]
    # Include exact immutable reference partitions and their original bindings.
    selected += [study/'clustering/baron-reference'/name for name in ['reference.json', 'partitions.npz']]
    selected += [study/'clustering/baron-protocol.json', study/'graph/baron-independent-check-final/checks.json']
    selected += [study/'clustering/independent-reference'/name for name in ['reference.json', 'partitions.npz', 'seed-7.json', 'seed-19.json', 'seed-41.json']]
    selected += [study/'clustering'/name for name in ['protocol.json', 'reference_clustering.py', 'reference-prepare-status.json', 'reference-prepare.log']]
    selected += [study/'graph/independent-check/checks.json', study/'graph/full/receipt.json', root/'full-reference-check/checks.json']
    selected += [study/'clustering/baron/execution.json', study/'clustering/baron/receipt.json']
    paths = [('study/'+str(path.relative_to(study)), path) for path in selected]
    for name in ['reference_clustering.py', 'check_reference_reuse.py', 'collect_clustering.py',
                 'requirements-clustering.txt', 'archive_clustering_reference.py', 'verify_clustering_reference.py']:
        paths.append(('tools/'+name, Path(__file__).with_name(name)))
    assert len({name for name, _ in paths}) == len(paths)
    args.output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, source in sorted(paths):
        assert source.is_file() and not source.is_symlink()
        raw = source.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        encoded = gzip.compress(raw, compresslevel=6, mtime=0)
        target = args.output/(name+'.gz')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(encoded)
        assert sha(source) == digest
        records.append(dict(sourcePath=name, sourceBytes=len(raw), sourceSHA256=digest,
                            storedPath=name+'.gz', storedBytes=len(encoded),
                            storedSHA256=hashlib.sha256(encoded).hexdigest(), gzipEncoded=True))
    embedded_sources = {str(source) for _, source in paths}
    external = {name: digest for name, digest in qualification['sourceIdentities'].items()
                if ('/clustering/baron/' in name or '/clustering-storage-adaptive/baron/' in name)
                and name not in embedded_sources}
    upstream = {}
    for folder in ['2026-09-10-graph', '2026-09-10-storage-access']:
        previous = read(Path(__file__).parent/'evidence'/folder/'manifest.json')
        for record in previous['records']:
            upstream[record['sourceSHA256']] = dict(archive=folder, sourcePath=record['sourcePath'])
        for name, identity in previous.get('fullSources', {}).items():
            upstream[identity['SHA256']] = dict(archive=folder, sourcePath=name)
    assert all(digest in upstream for digest in external.values()), 'Missing external restoration source'
    restoration = {name: dict(**upstream[digest], SHA256=digest) for name, digest in external.items()}
    manifest = dict(schemaVersion=1, records=records, externalOriginalBaronBundleFiles=external,
                    upstreamRestorationMap=restoration,
                    reusedCollectorPayloads={name: collection['remotePayloads'][name] for name in collection['reused']},
                    scope='Complete Baron direct and equivalence checks, seven expected rejections, remote collection and cross-cohort rejection. Original graph/PCA bundles remain in the earlier graph/storage archives and study paths, with exact external identities here. All three original full-HIRISA igraph partitions and their independent every-edge objective check are retained. Native HIRISA clustering/replay and biological qualification are not asserted. The initial overly specific test diagnostic failure is retained.')
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(status='archived', members=len(records),
                         storedBytes=sum(r['storedBytes'] for r in records),
                         manifestSHA256=sha(args.output/'manifest.json'))))


if __name__ == '__main__':
    main()
