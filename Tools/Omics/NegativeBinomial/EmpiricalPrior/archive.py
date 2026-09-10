#!/usr/bin/env python3
"""Package exact benchmark tables and logs; inventory large external sources.

Full native reports and source matrices remain external, with stored and logical
hashes. --verify checks every committed archive member without external data.
"""
import argparse, gzip, hashlib, json
from pathlib import Path

def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--verify', action='store_true')
a = p.parse_args()
if a.verify:
    manifest = json.loads((a.out/'manifest.json').read_text())
    for e in manifest['entries']:
        path = a.out/e['path']
        assert digest(path) == e['sha256'], path
        if 'logicalSHA256' in e:
            assert hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest() == e['logicalSHA256'], path
    actual = {str(x.relative_to(a.out)) for x in a.out.rglob('*') if x.is_file()}
    assert actual == {e['path'] for e in manifest['entries']} | {'manifest.json'}
    print(json.dumps(dict(status='all-archive-members-verified', files=len(manifest['entries']))))
    raise SystemExit(0)

assert a.root and not a.out.exists()
a.out.mkdir(parents=True)
entries = []; external = []
for path in sorted(a.root.rglob('*')):
    if not path.is_file(): continue
    rel = path.relative_to(a.root)
    # Keep complete source/export/native output outside Git. Record each exact
    # file instead of pretending that selected feature tables restore a report.
    large_external = (path.suffix in {'.h5ad', '.Rda'} or path.name == 'numivivo' or
                      rel.parts[0] == 'wire' and path.name not in {'export.json', 'session-info.txt'} or
                      rel.name == 'report.json.gz')
    # Author implementation is provenance, not redistributed project code.
    author_source = rel.parts[0] == 'source' and path.suffix == '.R'
    if large_external or author_source:
        item = dict(path=str(path), bytes=path.stat().st_size, sha256=digest(path))
        if rel.name == 'report.json.gz':
            item['logicalSHA256'] = hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest()
            item['remotePath'] = '/Users/n/numivivo-empirical-prior-20260910/'+str(rel)
        external.append(item)
        continue
    data = path.read_bytes()
    compress = path.suffix in {'.tsv', '.log', '.txt'} or path.suffix == '.json' and len(data) > 131072
    target = a.out/(str(rel)+'.gz' if compress else rel)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(gzip.compress(data, mtime=0) if compress else data)
    item = dict(path=str(target.relative_to(a.out)), bytes=target.stat().st_size, sha256=digest(target))
    if compress or path.suffix == '.gz':
        logical = data if compress else gzip.decompress(data)
        item['logicalSHA256'] = hashlib.sha256(logical).hexdigest()
    entries.append(item)
for rec in json.loads((a.root/'protocol.json').read_text())['records']:
    path=Path(rec['baselineReport']); assert digest(path)==rec['baselineReportSHA256']
    external.append(dict(path=str(path), sha256=digest(path), role='exact measured baseline report'))
(a.out/'manifest.json').write_text(json.dumps(dict(entries=entries, external=external,
    qualification='Exact tables and logs; full source matrices and native reports remain external. Empirical prior and conditional MAP checks; no posterior coverage or FDR calibration claim.',
    archiveScriptSHA256=digest(Path(__file__))), sort_keys=True, indent=2)+'\n')
print(json.dumps(dict(files=len(entries), bytes=sum(e['bytes'] for e in entries), externalFiles=len(external))))
