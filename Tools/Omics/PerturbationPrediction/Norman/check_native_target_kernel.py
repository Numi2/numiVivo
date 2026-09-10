#!/usr/bin/env python3
"""Compare frozen native target-kernel reports with pinned independent GO predictions."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from combinations import counts, load, metrics, sha, write, REFERENCE_SHA256
from go_transfer import independent_prediction, solve_weights
from run_native_target_kernel import ANNOTATIONS, COVERAGE

PRIMARY = 'cb7e7ec6c41f5322c1bb4cbcd91d4151e5c4d2676f7c4e71eb0d1e3fb0fb8ec9'
SHUFFLE = 'b45876c197a9ae550cd10eb6923b219be6488404e551992f803d6893109a36a2'
METHODS = ['noChange', 'meanSingleResponse', 'meanSupportedResponse', 'goRidgeFixed', 'shuffledGoRidgeFixed']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('native', 'reference', 'primary', 'shuffle', 'prepared', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    assert sha(a.reference) == REFERENCE_SHA256
    assert sha(a.primary/'predictions.npz') == PRIMARY and sha(a.shuffle/'predictions.npz') == SHUFFLE
    assert sha(a.prepared/'annotations.json') == ANNOTATIONS and sha(a.prepared/'coverage.json') == COVERAGE
    assert sha(a.primary/'folds.json') == json.loads((a.primary/'receipt.json').read_text())['foldsSHA256']
    native_receipt = json.loads((a.native/'receipt.json').read_text())
    assert native_receipt['status'] == 'passed'
    a.out.mkdir(parents=True, exist_ok=False)
    ref, frozen, shuffle = load(a.reference), load(a.primary/'predictions.npz'), load(a.shuffle/'predictions.npz')
    old_folds = json.loads((a.primary/'folds.json').read_text()); annotations = json.loads((a.prepared/'annotations.json').read_text())
    with gzip.open(a.native/'training.json.gz', 'rb') as stream: training_bytes = stream.read()
    import hashlib
    assert hashlib.sha256(training_bytes).hexdigest() == native_receipt['trainingSHA256']
    training = json.loads(training_bytes); targets = frozen['targetIDs'].tolist(); names = ref['conditions'].tolist()
    assert [t['id'] for t in training['selection']['targets']] == targets
    assert training['featureIDs'] == frozen['featureIDs'].tolist() == ref['feature_ids'].tolist()
    raw = np.vstack([ref['counts'][names.index(t)] for t in ['control', *targets]])
    matrix = training['matrix']
    for i, row in enumerate(raw):
        indices = np.flatnonzero(row); start, end = matrix['rowOffsets'][i:i+2]
        assert matrix['featureIndices'][start:end] == indices.tolist()
        assert matrix['counts'][start:end] == row[indices].tolist()
    cpm = counts(raw)/raw.sum(axis=1, keepdims=True)*1e6
    baseline = np.log1p(cpm[0]); del training, matrix, training_bytes
    comparisons, results, coverage, diagnostic_differences = [], [], [], []
    max_prediction = max_weights = max_independent = max_cpm = 0.0
    for i, target in enumerate(targets):
        fold = a.native/'folds'/target
        with gzip.open(fold/'report.json.gz', 'rb') as stream: report_bytes = stream.read()
        manifest = json.loads((fold/'manifest.json').read_text()); entries = {e['path']: e for e in manifest['files']}
        assert hashlib.sha256(report_bytes).hexdigest() == entries['prediction/report.json']['sha256']
        report = json.loads(report_bytes); q = report['queries'][0]
        assert report['featureIDs'] == ref['feature_ids'].tolist()
        train = [t for t in targets if t != target]
        assert manifest['trainingTargets'] == train == old_folds[i]['trainingTargets']
        supported = [t for t in train if t in annotations]
        assert report['supportedTrainingTargets'] == supported
        usable = [targets.index(t)+1 for t in supported]
        y = np.log1p(cpm[usable])-baseline
        ordinary_raw = dict(noChange=baseline, meanSingleResponse=baseline+(np.log1p(cpm[[j+1 for j,t in enumerate(targets) if t != target]])-baseline).mean(axis=0),
                            meanSupportedResponse=baseline+y.mean(axis=0))
        native = {b['method']: np.array(b['expression']) for b in report['baselines']}
        diagnostics = {b['method']: b['diagnostic'] for b in report['baselines']}
        coverage.append(dict(target=target, status=q['status']))
        if target in annotations:
            assert q['status'] == 'predicted' and supported == old_folds[i]['descriptorTrainingTargets']
            terms = [set(annotations[t]['terms']) for t in supported]; query = set(annotations[target]['terms'])
            kernel = np.array([[len(x & z)/len(x | z) for z in terms] for x in terms]); similarities = np.array([len(query & x)/len(query | x) for x in terms])
            weights, _ = solve_weights(kernel, similarities, 1)
            np.testing.assert_allclose(q['weights'], weights, atol=1e-12, rtol=1e-11)
            np.testing.assert_allclose(q['weights'], old_folds[i]['weights']['goRidgeFixed'], atol=1e-12, rtol=1e-11)
            max_weights = max(max_weights, float(np.max(abs(np.array(q['weights'])-weights))))
            for method, response, field, diagnostic in [('goRidgeFixed', y, 'expression', 'diagnostic'), ('shuffledGoRidgeFixed', np.roll(y,-1,axis=0), 'shuffledExpression', 'shuffledDiagnostic')]:
                ordinary_raw[method] = baseline+weights@response
                oracle = np.maximum(baseline+independent_prediction(kernel, similarities, response, 1), 0)
                native[method] = np.array(q[field]); diagnostics[method] = q[diagnostic]
                np.testing.assert_allclose(native[method], oracle, atol=1e-11, rtol=1e-11)
                max_independent = max(max_independent, float(np.max(abs(native[method]-oracle))))
        else:
            assert q['status'] in ('noData', 'unresolvedIdentity')
            assert all(q.get(k) is None for k in ('weights','expression','shuffledExpression','diagnostic','shuffledDiagnostic'))
        truth = np.log1p(counts(raw[i+1])/raw[i+1].sum()*1e6)-baseline
        for method, values in native.items():
            expected = shuffle['predictions'][i] if method == 'shuffledGoRidgeFixed' else frozen[method][i]
            if not np.all(np.isfinite(expected)): expected = np.maximum(ordinary_raw[method], 0)
            np.testing.assert_allclose(values, expected, atol=1e-11, rtol=1e-11)
            error = float(np.max(abs(values-expected))); max_prediction = max(max_prediction, error)
            d = diagnostics[method]
            assert d['renormalized'] is False
            cpmsum = float(np.expm1(values).sum()); cpm_error = abs(cpmsum-d['impliedCPMSum']); max_cpm = max(max_cpm, cpm_error)
            np.testing.assert_allclose(cpmsum, d['impliedCPMSum'], atol=1e-7, rtol=1e-12)
            external_clips = int(np.count_nonzero(ordinary_raw[method]<0))
            if external_clips != d['clippedFeatures']:
                ambiguous = int(np.count_nonzero(abs(ordinary_raw[method])<=1e-11))
                assert abs(external_clips-d['clippedFeatures']) <= ambiguous
                diagnostic_differences.append(dict(target=target, method=method, nativeClips=d['clippedFeatures'], independentClips=external_clips,
                                                   featuresWithinNumericalToleranceOfZero=ambiguous, maximumPredictionError=error))
            comparisons.append(dict(target=target, method=method, maximumPredictionError=error, maximumCPMSumError=cpm_error))
            for panel, ix in [('allGenes', np.arange(len(baseline))), ('trainingTop1000', frozen['trainingTop1000'][i])]:
                results.append(dict(target=target, supported=target in annotations, method=method, panel=panel, **metrics((values-baseline)[ix],truth[ix])))
    summary = []
    for cohort in ('allAvailable', 'matchedSupported'):
        for method in METHODS:
            for panel in ('allGenes','trainingTop1000'):
                rows = [r for r in results if r['method']==method and r['panel']==panel and (cohort=='allAvailable' or r['supported'])]
                zero = {r['target']:r['rmse'] for r in results if r['method']=='noChange' and r['panel']==panel}
                summary.append(dict(cohort=cohort, method=method, panel=panel, targets=len(rows), meanResponseRMSE=float(np.mean([r['rmse'] for r in rows])),
                                    worseThanNoChange=[r['target'] for r in rows if r['rmse']>zero[r['target']]]))
    write(a.out/'comparisons.json', comparisons); write(a.out/'results.json', results); write(a.out/'summary.json', summary)
    write(a.out/'diagnostic-differences.json', diagnostic_differences); write(a.out/'coverage.json', coverage)
    result = dict(status='passed', exactTrainingCounts=True, targets=105, predictedTargets=101, features=33694, comparedVectors=len(comparisons),
                  maximumFrozenPredictionError=max_prediction, maximumIndependentPredictionError=max_independent, maximumWeightError=max_weights,
                  maximumCPMSumError=max_cpm, clippingClassificationDifferences=len(diagnostic_differences), nativeReceiptSHA256=sha(a.native/'receipt.json'),
                  primaryPredictionSHA256=PRIMARY, shufflePredictionSHA256=SHUFFLE, referenceSHA256=REFERENCE_SHA256, checkerSHA256=sha(Path(__file__)))
    write(a.out/'checks.json', result); print(json.dumps(result))

if __name__ == '__main__':
    with threadpool_limits(limits=1): main()
