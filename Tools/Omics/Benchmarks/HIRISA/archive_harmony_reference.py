#!/usr/bin/env python3
"""Archive all complete HIRISA Harmony reference reports and response comparisons."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from archive_full_integration_reports import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    def read(path): return json.loads(path.read_text())
    complete = read(a.root/'complete.json'); protocol = read(a.root/'protocol.json')
    assert complete['status'] == 'three-complete-reference-fits-and-response-evaluations'
    assert complete['cells'] == protocol['cells'] == 1612594
    assert complete['seeds'] == protocol['seeds'] == [7, 19, 41]
    assert all(c['returnCode'] == 0 for c in complete['commands'])
    qualified = read(a.root/'qualification/checks.json')
    assert qualified['status'] == 'passed' and qualified['cells'] == 24673 and qualified['maximumAbsoluteError'] == 0
    external = []
    for seed in protocol['seeds']:
        case = a.root/f'seed-{seed}'; report = read(case/'report.json'); response = read(case/'response.json')
        assert report['status'] == 'measured' and report['cells'] == response['cells'] == 1612594
        assert report['bindings']['protocol']['SHA256'] == sha(a.root/'protocol.json')
        assert report['output'] == response['bindings']['scores']
        assert report['output']['SHA256'] == sha(case/'scores.bin')
        assert report['output']['bytes'] == (case/'scores.bin').stat().st_size
        assert report['runnerSHA256'] == sha(a.root/'run_harmony_integration.py')
        assert report['readerSHA256'] == sha(a.root/'check_integration_response.py')
        external.append(dict(sourcePath=f'seed-{seed}/scores.bin', **report['output']))
    a.output.mkdir(parents=True, exist_ok=False); records = []
    for source in sorted(a.root.rglob('*')):
        if not source.is_file() or '__pycache__' in source.parts or source.suffix not in ['.py', '.json', '.log']: continue
        assert not source.is_symlink()
        name = str(source.relative_to(a.root)); raw = source.read_bytes(); digest = hashlib.sha256(raw).hexdigest()
        encoded = gzip.compress(raw, compresslevel=6, mtime=0)
        target = a.output/(name+'.gz'); target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(encoded)
        assert sha(source) == digest
        records.append(dict(sourcePath=name, sourceBytes=len(raw), sourceSHA256=digest,
                            storedPath=name+'.gz', storedBytes=len(encoded),
                            storedSHA256=hashlib.sha256(encoded).hexdigest(), gzipEncoded=True))
    manifest = dict(schemaVersion=1, records=records, externalPayloads=external,
                    archiveScriptSHA256=sha(Path(__file__)),
                    scope='All three complete HIRISA Harmony2 reference fits and the original frozen response diagnostics. Full corrected matrices are preserved externally under rechecked exact identities. Same transductive PCA and experimental matching; four insensitive classification controls remain insufficient, and one native seed does not establish native multiseed robustness or broad biological validation.')
    (a.output/'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(members=len(records), storedBytes=sum(r['storedBytes'] for r in records),
                         manifestSHA256=sha(a.output/'manifest.json'))))


if __name__ == '__main__':
    main()
