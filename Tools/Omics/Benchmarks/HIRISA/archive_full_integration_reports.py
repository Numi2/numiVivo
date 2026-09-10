#!/usr/bin/env python3
"""Freeze complete HIRISA integration reports, with large payloads bound externally."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'response', 'output']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    def read(path): return json.loads(path.read_text())
    publication = read(a.root/'full-output-freeze.json')
    numerical = read(a.root/'independent-check/checks.json')
    link = read(a.response/'native-evaluation-link.json')
    response = read(a.response/link['resultFile'])
    protocol = read(a.response/'protocol.json')
    assert publication['status'] == 'published'
    assert read(a.root/'full-publish-status.json')['returnCode'] == 0
    assert numerical['status'] == 'passed' and numerical['cells'] == protocol['cells'] == 1612594
    assert read(a.root/'independent-check-status.json')['returnCode'] == 0
    assert sha(a.response/link['resultFile']) == link['resultSHA256']
    assert link['nativePublicationSHA256'] == sha(a.root/'full-output-freeze.json')
    assert response['bindings']['scores'] == publication['files']['scores.bin']
    assert response['bindings']['metadata'] == publication['files']['metadata.json']
    assert response['bindings']['protocol']['SHA256'] == sha(a.response/'protocol.json')
    assert response['cells'] == 1612594 and len(response['folds']) == 79
    # Report measurements honestly; a preservation failure is still an outcome
    # to archive, not a reason to discard the experimental record.
    external = []
    for name, expected in publication['files'].items():
        source = a.root/'native'/name
        assert source.stat().st_size == expected['bytes'] and sha(source) == expected['SHA256']
        if name not in ['report.json', 'plan.json', 'receipt.json']:
            external.append(dict(sourcePath='integration-full/native/'+name, **expected))
    selected = ['protocol.json', 'plan.json', 'pca-plan.json', 'pca-fit-status.json', 'pca-exact.json',
                'full-publish-start.json', 'full-publish-status.json', 'full-publish.log',
                'full-output-freeze.json', 'full-publish-admission-refused.json',
                'reference-not-started.json', 'resume-start.json', 'reference-resume-start.json',
                'reference-execution-freeze.json', 'independent-check-start.json',
                'independent-check-status.json', 'independent-check.log', 'independent-check/checks.json',
                'native/report.json', 'native/plan.json', 'native/receipt.json']
    native_replay = 'pending'
    if (a.root/'full-verify-status.json').exists():
        status = read(a.root/'full-verify-status.json')
        native_replay = 'passed' if status['returnCode'] == 0 else 'failed'
        selected += ['full-verify-start.json', 'full-verify-status.json', 'full-verify.log']
        if (a.root/'native-complete.json').exists(): selected.append('native-complete.json')
    paths = [('integration-full/'+name, a.root/name) for name in selected]
    paths += [('integration-response-v2/'+path.name, path) for path in sorted(a.response.iterdir())
              if path.is_file() and path.suffix in ['.json', '.py', '.log']]
    a.output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, source in paths:
        assert source.is_file() and not source.is_symlink()
        raw = source.read_bytes(); before = hashlib.sha256(raw).hexdigest()
        encoded = gzip.compress(raw, compresslevel=6, mtime=0)
        target = a.output/(name+'.gz'); target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(encoded)
        assert sha(source) == before
        records.append(dict(sourcePath=name, sourceBytes=len(raw), sourceSHA256=before,
                            storedPath=name+'.gz', storedBytes=len(encoded),
                            storedSHA256=hashlib.sha256(encoded).hexdigest(), gzipEncoded=True))
    manifest = dict(schemaVersion=1, records=records, externalPayloads=external,
                    nativeReplay=native_replay, archiveScriptSHA256=sha(Path(__file__)),
                    allResponseDiagnosticGatesPassed=response['allResponseDiagnosticGatesPassed'],
                    classificationControlSensitiveContrasts=12, classificationControlInsufficientContrasts=4,
                    scope='Complete original HIRISA native integration reports and independent numerical/response diagnostics. All published payloads rehashed; large matrices and metadata remain preserved externally under exact identities, not embedded in this report archive. Four insensitive classification controls and prior missing-design launch failure retained. Coarse response diagnostics do not establish complete biological preservation or independent-method superiority.')
    (a.output/'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(members=len(records), storedBytes=sum(r['storedBytes'] for r in records),
                         nativeReplay=native_replay, manifestSHA256=sha(a.output/'manifest.json'))))


if __name__ == '__main__':
    main()
