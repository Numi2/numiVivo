#!/usr/bin/env python3
"""Archive source/design/interop evidence, leaving complete source payloads external."""
import argparse
import gzip
import json
from pathlib import Path

from acquire import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--allow-native-pending', action='store_true')
    args = parser.parse_args()
    root, output = args.root, args.output
    assert json.loads((root / 'source-audit/summary.json').read_text())['complete']
    assert json.loads((root / 'anndata-verification.json').read_text())['status'] == 'passed'
    native_complete = (root / 'native-baseline-status.json').is_file()
    assert native_complete or args.allow_native_pending
    native_status_path = root / ('native-baseline-status.json' if native_complete else 'native-pending.json')
    native_status = json.loads(native_status_path.read_text())
    if not native_complete:
        assert native_status['status'] == 'transfer-running' and not native_status['nativeExecutionQualified']
    output.mkdir(parents=True, exist_ok=False)
    names = ['series.soft.gz', 'filelist.txt', 'design.json', 'source-receipt.json',
             'prepared-receipt.json', 'stream-plan.json', 'anndata-verification.json',
             'acquire.log', 'audit-initial.log', 'audit-second.log', 'audit-full.log',
             'audit-full-status.json', 'prepare-full.log', 'prepare-full-status.json',
             'anndata-verification.log', 'anndata-verification-status.json',
             'baseline-capacity-rejection.json']
    names += ['native-baseline.log', 'native-baseline-status.json', 'baseline-transport.log'] if native_complete else ['native-pending.json']
    names += [str(p.relative_to(root)) for p in sorted((root / 'source-audit').iterdir())
              if p.suffix in {'.json', '.npz'}]
    records = []
    for name in names:
        path = root / name
        data = path.read_bytes()
        suffix = '' if path.suffix in {'.gz', '.npz'} else '.gz'
        destination = output / (name + suffix)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(gzip.compress(data, mtime=0) if suffix else data)
        records.append({'sourcePath': name, 'storedPath': str(destination.relative_to(output)),
                        'sourceBytes': len(data), 'sourceSHA256': digest(path),
                        'storedBytes': destination.stat().st_size, 'storedSHA256': digest(destination),
                        'gzipEncoded': bool(suffix)})
    manifest = {'schemaVersion': 1, 'records': records,
                'sourceFiles': json.loads((root / 'source-receipt.json').read_text())['files'],
                'largeSourcePayloadsExternal': True,
                'nativeStatus': native_status}
    (output / 'manifest.json').write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'files': len(records), 'bytes': sum(r['storedBytes'] for r in records),
                      'manifestSHA256': digest(output / 'manifest.json')}))

if __name__ == '__main__':
    main()
