#!/usr/bin/env python3
"""Archive compact global-scale evidence; bind large matrices by exact hashes."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path)
parser.add_argument('--out', type=Path, required=True)
parser.add_argument('--verify', action='store_true')
args = parser.parse_args()


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


if args.verify:
    manifest = json.loads((args.out / 'manifest.json').read_text())
    for entry in manifest['entries']:
        raw = (args.out / entry['path']).read_bytes()
        assert sha(raw) == entry['sha256'], entry['path']
        if 'logicalSHA256' in entry:
            assert sha(gzip.decompress(raw)) == entry['logicalSHA256'], entry['path']
    assert {str(p.relative_to(args.out)) for p in args.out.rglob('*') if p.is_file()} == {
        e['path'] for e in manifest['entries']} | {'manifest.json'}
    print(json.dumps(dict(status='all-archive-members-verified', files=len(manifest['entries']))))
    raise SystemExit(0)

assert args.root and not args.out.exists()
summary = json.loads((args.root / 'summary.json').read_text())
assert summary['status'] == 'passed' and summary['arms'] == 58
args.out.mkdir(parents=True)
entries, external = [], []
large_names = {'native.json.gz', 'reference.json.gz', 'independent-input.json.gz', 'independent.json.gz'}
for path in sorted(args.root.rglob('*')):
    if not path.is_file():
        continue
    assert not path.is_symlink(), path
    relative = path.relative_to(args.root)
    raw = path.read_bytes()
    if path.name in large_names:
        external.append(dict(path=str(path), bytes=len(raw), sha256=sha(raw),
                             logicalSHA256=sha(gzip.decompress(raw)),
                             role='Full scale-study input, native or independent reference matrices'))
        continue
    compress = path.suffix in {'.log', '.txt', '.swift', '.R'} or (
        path.suffix == '.json' and len(raw) > 32768)
    data = gzip.compress(raw, mtime=0) if compress else raw
    destination = args.out / (str(relative) + '.gz' if compress else str(relative))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(data)
    entry = dict(path=str(destination.relative_to(args.out)), bytes=len(data), sha256=sha(data))
    if compress or path.suffix == '.gz':
        entry['logicalSHA256'] = sha(raw if compress else gzip.decompress(raw))
    entries.append(entry)
execution = json.loads((args.root / 'execution.json').read_text())
external.extend(execution['externalExecutables'])
external.extend(execution['priorInputs'])
manifest = dict(entries=entries, external=external, protocolSHA256=execution['protocolSHA256'],
                qualification='58 complete native global-scale/refit arms; independent native-input LOWESS and fixed-scale GLMs; recomputed scores; controlled LOWESS, Swift tests and retained failed checker/compiler attempts. Large matrices remain external. No QL prior, moderation, cohort inference or calibration claim.')
(args.out / 'manifest.json').write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
print(json.dumps(dict(files=len(entries), storedBytes=sum(e['bytes'] for e in entries),
                      externalFiles=len(external))))
