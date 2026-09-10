#!/usr/bin/env python3
"""Verify complete global-scale evidence and independently recompute NB scores."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read(path, expected=None):
    raw = path.read_bytes()
    if expected is not None:
        assert sha(raw) == expected, path
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def scores(counts, means, design, dispersion):
    denominator = 1 + dispersion[:, None] * means
    return np.max(np.abs(((counts - means) / denominator) @ design) /
                  np.sqrt((means / denominator) @ (design * design)), axis=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--ql-root', type=Path, required=True)
    args = parser.parse_args()
    root, source = args.root, Path(__file__).parent
    protocol_sha = sha((source / 'PROTOCOL.md').read_bytes())
    protocol = read(root / 'native/protocol.json')
    assert protocol['protocolSHA256'] == protocol_sha
    expected = {(c['id'], m) for c in protocol['cases']
                for m in ('nativeTrend-adjusted', 'edgeRTrend-adjusted')}
    runs = read(root / 'native/complete.json')
    independent = read(root / 'native/independent-complete.json')
    assert len(expected) == len(runs) == len(independent) == 58
    assert {(r['case'], r['method']) for r in runs} == expected
    assert {(r['case'], r['method']) for r in independent} == expected
    independent = {(r['case'], r['method']): r for r in independent}
    lowess = read(root / 'lowess/checks.json')
    assert lowess['status'] == 'passed' and len(lowess['checks']) == 25
    assert lowess['protocolSHA256'] == protocol_sha
    results, failures = [], []
    for run in runs:
        case, method = run['case'], run['method']
        dest = root / 'native' / case / method
        check = independent[(case, method)]
        assert read(dest / 'run.json') == run
        assert read(dest / 'independent-receipt.json') == check
        assert run['binarySHA256'] == protocol['binarySHA256']
        assert run['status'] == 'passed-native-stage' and run['exitCode'] == 0
        assert check['status'] == 'passed' and check['exitCode'] == 0
        assert check['nativeSHA256'] == run['outputSHA256']
        assert check['scriptSHA256'] == sha((source / 'check_reference.R').read_bytes())
        data = read(args.ql_root / case / 'input.json.gz', run['inputSHA256'])
        ref_receipt = read(root / 'reference' / case / 'reference-receipt.json')
        assert ref_receipt['inputSHA256'] == run['inputSHA256']
        assert ref_receipt['outputSHA256'] == run['referenceSHA256']
        assert ref_receipt['protocolSHA256'] == protocol_sha
        assert ref_receipt['scriptSHA256'] == sha((source / 'reference.R').read_bytes())
        reference = read(root / 'reference' / case / 'reference.json.gz',
                         run['referenceSHA256'])['results'][method]
        prior = read(args.ql_root / case / 'reference.json.gz',
                     ref_receipt['priorSHA256'])['results'][method]
        assert reference['priorScaleDifference'] <= 1e-10
        assert reference['priorMeansDifference'] <= 1e-8
        compact = read(dest / 'independent-input.json.gz', check['requestSHA256'])
        comparison = read(dest / 'independent.json.gz', check['outputSHA256'])
        native = read(dest / 'native.json.gz', run['outputSHA256'])['fit']
        assert native['completed'] and not native['failures']
        assert len(native['updates']) == 2
        assert native['updates'][0]['inputScale'] == 1
        assert native['updates'][1]['inputScale'] == native['updates'][0]['outputScale']
        scale = native['averageQuasiDispersion']
        assert scale == native['updates'][1]['outputScale'] and scale >= 1
        assert compact['averageQuasiDispersion'] == scale
        assert compact['trendDispersions'] == reference['trendDispersions']
        assert compact['abundanceCovariates'] == reference['abundanceCovariates']
        y, x = np.asarray(data['counts'], dtype=float), np.asarray(data['design'])
        phi = np.asarray(compact['trendDispersions'])
        assert len(data['featureIndices']) == run['genes'] == len(native['adjustedResiduals'])
        score_values, score_errors = [], []
        for label, fits, dispersion in (
                ('initial', native['initialFits'], phi),
                ('final', native['refittedFits'], phi / scale)):
            assert [f['means'] for f in fits] == compact[label + 'Means']
            assert [f['effect'] for f in fits] == compact[label + 'Effects']
            assert all(f['converged'] and not f['positiveCountDesignRankDeficient'] for f in fits)
            assert np.array_equal([f['dispersion'] for f in fits], dispersion)
            score = scores(y, np.asarray(compact[label + 'Means']), x, dispersion)
            discrepancy = np.abs(score - np.asarray([f['maximumScaledScore'] for f in fits]))
            # 1e-10 absolute allowance covers independent floating-point arithmetic,
            # not a new native solver stopping threshold.
            if not np.all(np.isfinite(score)) or np.max(score) > 1e-7 + 1e-10:
                failures.append(dict(case=case, method=method, stage=label,
                                     maximumScore=float(np.max(score))))
            if np.max(discrepancy) > 1e-10:
                failures.append(dict(case=case, method=method, stage=label,
                                     maximumReportedScoreDiscrepancy=float(np.max(discrepancy))))
            score_values.append(score)
            score_errors.append(float(np.max(discrepancy)))
        if run['maximumScaledScore'] > 1e-7:
            failures.append(dict(case=case, method=method, stage='native-reported-score'))
        table = comparison['table']
        assert len(table) == run['genes']
        residual_errors = {'deviance': 0.0, 'degreesOfFreedom': 0.0}
        for index, (feature, row, residual, old) in enumerate(zip(
                data['featureIndices'], table, native['adjustedResiduals'], prior['table'], strict=True)):
            assert row['featureIndex'] == feature == old['featureIndex']
            row['initialNativeScore'] = float(score_values[0][index])
            row['finalNativeScore'] = float(score_values[1][index])
            row['nativeResidualDeviance'] = residual['deviance']
            row['nativeResidualDF'] = residual['degreesOfFreedom']
            row['nativeQuasiDispersion'] = residual.get('quasiDispersion')
            row['originalReferenceResidualDeviance'] = old['residualDeviance']
            row['originalReferenceResidualDF'] = old['residualDF']
            for key, old_key in [('deviance', 'residualDeviance'), ('degreesOfFreedom', 'residualDF')]:
                error = abs(residual[key] - old[old_key]) / max(1, abs(old[old_key]))
                residual_errors[key] = max(residual_errors[key], error)
        smooth_checks = comparison['smootherChecks']
        assert len(smooth_checks) == 2
        for update, sm, request in zip(native['updates'], smooth_checks, compact['updates'], strict=True):
            assert update['smoother']['emptyWeightAnchors'] == 0
            assert update['eligibleIndices'] == request['eligibleIndices']
            assert update['quasiDispersions'] == request['quasiDispersions']
            assert update['smoother']['fitted'] == request['smoother']
            assert sorted(update['eligibleIndices'] + update['omittedIndices']) == list(range(run['genes']))
            predicted_scale = max(1, float(np.quantile(update['smoother']['fitted'], .9))) ** 4
            assert abs(predicted_scale / update['outputScale'] - 1) <= 1e-12
            assert sm['maximumRelativeError'] <= 2e-7
        max_fit_error = max(row[k] for row in table for k in (
            'initialMeanError', 'finalMeanError', 'initialEffectError', 'finalEffectError'))
        assert max_fit_error == check['maximumFitRelativeError'] <= 2e-5
        assert not any(row['initialReferenceFailed'] or row['finalReferenceFailed'] for row in table)
        compact_result = dict(table=table, smootherChecks=smooth_checks,
                              nativeUpdates=native['updates'])
        packed = gzip.compress(json.dumps(compact_result, separators=(',', ':'),
                                          allow_nan=False).encode(), mtime=0)
        (dest / 'checks.json.gz').write_bytes(packed)
        result = dict(case=case, method=method, genes=run['genes'],
                      nativeSHA256=run['outputSHA256'], independentSHA256=check['outputSHA256'],
                      checksSHA256=sha(packed), scale=scale, referenceScale=run['referenceScale'],
                      scaleRelativeDifference=run['scaleRelativeDifference'],
                      maximumFitRelativeError=max_fit_error,
                      maximumSmootherRelativeError=check['maximumSmootherRelativeError'],
                      maximumNativeScore=float(max(np.max(v) for v in score_values)),
                      maximumScoreDiscrepancy=max(score_errors),
                      maximumResidualRelativeDifference=residual_errors,
                      omittedIndices=run['omittedIndices'], momentEvaluations=run['momentEvaluations'],
                      wallSecondsIncludingTransport=run['seconds'])
        results.append(result)
        print('verified', case, method, flush=True)
    summary = dict(status='passed' if not failures else 'failed', arms=len(results),
                   geneArmFits=sum(r['genes'] for r in results), lowessCases=25,
                   controlledLowessMaximumRelativeError=lowess['maximumRelativeError'],
                   momentEvaluations=sum(r['momentEvaluations'] for r in results),
                   omittedGeneUpdates=sum(len(v) for r in results for v in r['omittedIndices']),
                   failures=failures, results=results)
    for key in ['maximumFitRelativeError', 'maximumSmootherRelativeError', 'maximumNativeScore',
                'maximumScoreDiscrepancy', 'scaleRelativeDifference']:
        summary[key] = max(r[key] for r in results)
    (root / 'summary.json').write_text(json.dumps(summary, sort_keys=True, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'results'}, indent=2))
    return 0 if not failures else 1


if __name__ == '__main__':
    raise SystemExit(main())
