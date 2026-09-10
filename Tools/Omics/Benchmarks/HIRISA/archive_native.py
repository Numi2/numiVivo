#!/usr/bin/env python3
"""Archive complete native ingestion gates, independently of ongoing DE/prediction."""
import argparse
import gzip
import json
import shutil
from pathlib import Path
from acquire import digest


def read(path):
    return json.loads(path.read_text())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args();root = a.root;out = a.output
    independent = read(root / 'native-independent-verification.json')
    assert independent['status'] == 'passed' and independent['nativeRawIngestionLinked']
    identity = read(root / 'native-report-identity.json')
    jobs = read(root / 'native-jobs-at-archival.json')
    assert identity['status'] == 'passed' and identity['sameReportBytes']
    assert identity['reportSHA256'] == independent['reportSHA256'] == digest(root / 'native-full/report.json')
    assert read(root / 'native-release-publish-status.json')['returnCode'] == 0
    assert read(root / 'native-release-verify-status.json')['returnCode'] == 0
    assert read(root / 'native-baseline-status.json')['exitCode'] == 65
    assert read(root / 'reader-checks/checks.json')['status'] == 'passed'
    assert read(root / 'reader-release-checks/checks.json')['status'] == 'passed'
    assert read(root / 'prediction-folds.json')['foldCount'] == 79
    assert '4 tests in 1 suite passed' in (root / 'snapshot-tests.log').read_text()
    names = ['native-baseline.log', 'native-baseline-status.json', 'baseline-capacity-rejection.json',
             'native-full-publish.log', 'native-full-publish-status.json', 'native-full-environment.json',
             'native-release-publish.log', 'native-release-publish-status.json',
             'native-release-verify.log', 'native-release-verify-status.json', 'native-release-environment.json',
             'native-report-identity.json', 'capture_native_identity.py', 'native-jobs-at-archival.json',
             'native-independent-verification.json', 'native-independent-verification.log',
             'prepare-native-qc.log', 'prediction-folds.json', 'freeze-prediction.log',
             'freeze_release_runtime.py', 'freeze-release-runtime.log', 'run_native_full.py', 'run_native_release.py',
             'remote_product.py', 'remote_product_release.py', 'reader-checks.log', 'reader-release-checks.log',
             'cli-build.log', 'cli-clone-build.log', 'cli-release-build.log', 'snapshot-tests.log',
             'native-full-inflight.sample.txt', 'native-full-inflight-sample.log',
             'prepared-receipt.json', 'source-receipt.json', 'design.json', 'stream-plan.json',
             'anndata-verification.json', 'inference-inputs/receipt.json', 'inference-inputs/counts.npz',
             'prepare-inference-attempt-1.py', 'prepare-inference-attempt-1.log', 'prepare-inference.log',
             'inference-transfer-attempt-1.json', 'model-check-available-3.log', 'archive-native-attempt-1.log']
    debug_replay = root / 'native-full-verify-status.json'
    if debug_replay.exists():
        names += ['native-full-verify-status.json', 'native-full-verify.log']
    for folder in ['native-full', 'native-release', 'native-qc-reference', 'reader-checks',
                   'reader-release-checks', 'inference-inputs/cohorts', 'inference-inputs/requests',
                   'inference-inputs/reference-inputs']:
        names += [str(path.relative_to(root)) for path in sorted((root / folder).rglob('*'))
                  if path.is_file() and not (folder in ['native-full', 'native-release'] and path.name == 'original.h5ad')]
    names = sorted(set(names))
    missing = [name for name in names if not (root / name).is_file()]
    assert not missing, missing
    out.mkdir(parents=True, exist_ok=False);records = []
    for name in names:
        source = root / name
        assert not source.is_symlink()
        encoded = source.suffix not in ['.gz', '.npz', '.h5ad']
        target = out / (name + ('.gz' if encoded else ''))
        target.parent.mkdir(parents=True, exist_ok=True)
        before = digest(source)
        with source.open('rb') as src, target.open('xb') as dst:
            if encoded:
                with gzip.GzipFile(fileobj=dst, mode='wb', filename='', mtime=0) as zipped:
                    shutil.copyfileobj(src, zipped, 1_048_576)
            else:
                shutil.copyfileobj(src, dst, 1_048_576)
        assert digest(source) == before, ('source changed during archival', name)
        records.append({'sourcePath': name, 'storedPath': str(target.relative_to(out)),
                        'sourceBytes': source.stat().st_size, 'sourceSHA256': before,
                        'storedBytes': target.stat().st_size, 'storedSHA256': digest(target),
                        'gzipEncoded': encoded})
    manifest = {'schemaVersion': 1, 'records': records, 'sourceSHA256': independent['sourceSHA256'],
                'nativeIndependentVerification': independent, 'crossBuildReportIdentity': identity,
                'releaseNativeReplay': read(root / 'native-release-verify-status.json'),
                'debugNativeReplay': read(debug_replay) if debug_replay.exists() else None,
                'debugReplayTerminalAtArchival': debug_replay.exists(),
                'observedJobsAtCapture': jobs,
                'largeSourcePayloadsExternal': True, 'nativeExecutablesRetainedExternally': True,
                'pairedDEAndPredictionBenchmarkComplete': False,
                'scope': 'Complete-source ingestion, independent cell/QC/count linkage and release replay. Original 131 sources, full H5AD and executables remain external. DE output families are a separate ongoing archive; prediction folds are frozen and unfitted.'}
    (out / 'manifest.json').write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'members': len(records), 'storedBytes': sum(x['storedBytes'] for x in records),
                      'manifestSHA256': digest(out / 'manifest.json')}))


if __name__ == '__main__':
    main()
