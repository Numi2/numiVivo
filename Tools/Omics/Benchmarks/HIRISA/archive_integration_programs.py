#!/usr/bin/env python3
"""Archive complete measured HIRISA program diagnostics, including failed gates."""
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


def read(path):
    return json.loads(path.read_text())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    root = a.study/'integration-programs'
    complete = read(root/'complete.json')
    assert complete['status'] == 'seven-complete-program-evaluations'
    assert len(complete['commands']) == 7
    assert all(v['returnCode'] == 0 for v in complete['commands'])
    oracles = read(root/'dense-oracles-complete.json')
    assert oracles['status'] == 'passed' and len(oracles['commands']) == 5
    assert all(v['returnCode'] == 0 for v in oracles['commands'])
    for name, digest in read(root/'evaluation-execution-freeze.json')['files'].items():
        assert sha(root/name) == digest, name
    external = {}
    for name in ['baseline', 'native', 'harmony-7', 'harmony-19', 'harmony-41']:
        result = read(root/(name+'.json'))
        oracle = read(root/(name+'-dense-oracle.json'))
        assert oracle['status'] == 'passed' and oracle['folds'] == 79
        assert oracle['resultSHA256'] == sha(root/(name+'.json'))
        assert result['cells'] == 1612594 and result['matchedOriginalCells'] == 1275710
        external[name] = {k: result['bindings'][k] for k in ['source', 'metadata', 'scores']}
    # These bindings identify existing external source/latent bundles. The evaluator
    # checked them before and after execution; this archiver does not claim to
    # re-read external matrices or to embed them in the report archive.
    paths = [('programs/'+str(x.relative_to(root)), x) for x in root.rglob('*')
             if x.is_file() and '__pycache__' not in x.parts
             and x.suffix in {'.json', '.log', '.py', '.npy'}]
    for rel in ['integration-full/full-output-freeze.json', 'integration-full/native/receipt.json',
                'native-qc-reference/receipt.json']:
        paths.append(('study/'+rel, a.study/rel))
    for seed in [7, 19, 41]:
        rel = f'integration-harmony/seed-{seed}/report.json'
        paths.append(('study/'+rel, a.study/rel))
    for name in ['prepare_program_reference.py', 'check_integration_programs.py',
                 'check_integration_response.py', 'test_integration_programs.py',
                 'check_program_dense_oracle.py', 'archive_integration_programs.py']:
        paths.append(('tools/'+name, Path(__file__).with_name(name)))
    assert len({n for n, _ in paths}) == len(paths)
    for _, source in paths:
        assert source.is_file() and not source.is_symlink(), source
    a.output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, source in sorted(paths):
        assert source.is_file() and not source.is_symlink(), source
        raw = source.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        encoded = gzip.compress(raw, compresslevel=6, mtime=0)
        target = a.output/(name+'.gz')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(encoded)
        assert sha(source) == digest
        records.append(dict(sourcePath=name, sourceBytes=len(raw), sourceSHA256=digest,
                            storedPath=name+'.gz', storedBytes=len(encoded),
                            storedSHA256=hashlib.sha256(encoded).hexdigest(), gzipEncoded=True))
    manifest = dict(schemaVersion=1, records=records, externalBindings=external,
                    scope='All original cells and program arrays; seven full diagnostics and five independent weighted-SVD checks. Failed and insensitive gates retained. External original H5AD, metadata and latent matrices remain in study storage with exact identities. This is an independent RNA reference and integration diagnostic, not newly executed native million-cell program scoring or prospective phenotype prediction.')
    (a.output/'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(status='archived', members=len(records),
                         storedBytes=sum(v['storedBytes'] for v in records),
                         manifestSHA256=sha(a.output/'manifest.json'))))


if __name__ == '__main__':
    main()
