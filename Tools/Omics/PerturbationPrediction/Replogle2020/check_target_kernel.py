#!/usr/bin/env python3
"""Independently check every frozen native prediction, then score all five groups."""
import argparse
import gzip
import hashlib
import importlib.metadata
import json
import sys
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'Norman'))
from combinations import metrics
from go_transfer import independent_prediction, solve_weights

from prepare_training import dense_row, sha, write

REFERENCE = '97415859bbaa2c52c9244e45a3bd1d4c864abfad777eed36c08bf88773f61442'
METHODS = ['noChange', 'meanSingleResponse', 'meanSupportedResponse', 'goRidgeFixed', 'shuffledGoRidgeFixed']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--descriptors', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root, native, descriptors, out = args.root, args.native, args.descriptors, args.out
    frozen = json.loads((native/'prediction-freeze.json').read_text())
    assert frozen['status'] == 'all-native-predictions-frozen' and frozen['folds'] == 150
    assert sha(native/'input-freeze.json') == frozen['inputFreezeSHA256']
    assert sha(native/'folds.json') == frozen['foldManifestSHA256']
    inputs = json.loads((native/'input-freeze.json').read_text())
    assert sha(descriptors/'receipt.json') == inputs['descriptorReceiptSHA256']
    assert sha(descriptors/'annotations.json') == inputs['annotationsSHA256']
    assert sha(root/'selected-reference.npz') == REFERENCE
    assert frozen['allNativeReplaysPassed'] and frozen['firstRepeatedPredictionExact']
    out.mkdir(parents=True, exist_ok=False)
    reference = np.load(root/'selected-reference.npz')
    inventory = json.loads((root/'prepared-complete/inventory.json').read_text())
    source_rows = {name: i for i, name in enumerate(inventory['groups'])}
    annotations = json.loads((descriptors/'annotations.json').read_text())
    manifests = {(x['gemgroup'], x['target']): x for x in json.loads((native/'folds.json').read_text())}
    count_freeze = json.loads((root/'training-preparation/count-folds.json').read_text())
    raw_by_context, features_by_context, targets_by_context = {}, {}, {}
    for context in count_freeze['contexts']:
        gem = context['gemgroup']
        raw_bytes = gzip.decompress((native/(gem+'-training.json.gz')).read_bytes())
        assert hashlib.sha256(raw_bytes).hexdigest() == context['trainingSHA256']
        training = json.loads(raw_bytes)
        targets = [x['id'] for x in training['selection']['targets']]
        control = sum((reference['counts'][source_rows[gem+'|'+g]].astype(np.uint64)
                       for g in ['sgNegCtrl2', 'sgNegCtrl3']), np.zeros(33694, dtype=np.uint64))
        raw = np.vstack([control, *[reference['counts'][source_rows[gem+'|'+g]] for g in targets]])
        assert raw.shape == (31, 33694) and np.all(raw.sum(axis=1) > 0)
        for row in range(31):
            np.testing.assert_array_equal(dense_row(training['matrix'], row), raw[row])
        raw_by_context[gem] = raw.astype(np.float64)
        features_by_context[gem] = training['featureIDs']
        targets_by_context[gem] = targets
    comparisons, scores, diagnostic_differences, panels = [], [], [], []
    for item in inputs['folds']:
        gem, target = item['gemgroup'], item['target']
        fold = native/'folds'/gem/target
        assert sha(fold/'plan.json') == item['planSHA256']
        assert sha(fold/'query-input.json') == item['querySHA256']
        manifest = manifests[gem, target]
        report_bytes = gzip.decompress((fold/'report.json.gz').read_bytes())
        assert hashlib.sha256(report_bytes).hexdigest() == next(x['SHA256'] for x in manifest['files'] if x['path']=='prediction/report.json')
        report = json.loads(report_bytes)
        assert report['featureIDs'] == features_by_context[gem]
        targets = targets_by_context[gem]
        held = targets.index(target)
        keep = [i for i, t in enumerate(targets) if t != target]
        train = [targets[i] for i in keep]
        assert manifest['trainingTargets'] == train and item['keep'] == keep
        raw = raw_by_context[gem]
        cpm = raw/raw.sum(axis=1, keepdims=True)*1e6
        baseline = np.log1p(cpm[0])
        delta = np.log1p(cpm[np.asarray(keep)+1])-baseline
        supported = [t for t in train if t in annotations]
        assert report['supportedTrainingTargets'] == supported
        positions = [train.index(t) for t in supported]
        response = delta[positions]
        expected_raw = dict(noChange=baseline, meanSingleResponse=baseline+delta.mean(axis=0),
                            meanSupportedResponse=baseline+response.mean(axis=0))
        values = {b['method']: np.asarray(b['expression']) for b in report['baselines']}
        diagnostics = {b['method']: b['diagnostic'] for b in report['baselines']}
        query = report['queries'][0]
        maximum_weight_error = 0.0
        if target in annotations:
            assert query['status'] == 'predicted'
            terms = [set(annotations[t]['terms']) for t in supported]
            query_terms = set(annotations[target]['terms'])
            kernel = np.asarray([[len(x & y)/len(x | y) for y in terms] for x in terms])
            similarities = np.asarray([len(query_terms & t)/len(query_terms | t) for t in terms])
            weights, _ = solve_weights(kernel, similarities, 1)
            np.testing.assert_allclose(query['weights'], weights, atol=1e-12, rtol=1e-11)
            maximum_weight_error = float(np.max(abs(np.asarray(query['weights'])-weights)))
            for method, y, field, diagnostic in [
                ('goRidgeFixed', response, 'expression', 'diagnostic'),
                ('shuffledGoRidgeFixed', np.roll(response, -1, axis=0), 'shuffledExpression', 'shuffledDiagnostic')]:
                expected_raw[method] = baseline+weights@y
                oracle = np.maximum(baseline+independent_prediction(kernel, similarities, y, 1), 0)
                values[method] = np.asarray(query[field])
                diagnostics[method] = query[diagnostic]
                np.testing.assert_allclose(values[method], oracle, atol=1e-11, rtol=1e-11)
        else:
            assert query['status'] in ['noData', 'unresolvedIdentity']
        mean_cpm = (cpm[np.asarray(keep)+1].sum(axis=0)+cpm[0])/(len(keep)+1)
        top = np.asarray(sorted(range(33694), key=lambda i: (-mean_cpm[i], features_by_context[gem][i]))[:1000])
        panels.append(dict(gemgroup=gem, target=target, trainingTop1000=top.tolist()))
        # Held expression enters metrics only after the complete native freeze.
        # Use the same verified native control vector for both response differences.
        # This makes the declared no-change response exactly zero instead of
        # turning cross-library log1p rounding into a spurious sign/correlation.
        response_baseline = values['noChange']
        truth = np.log1p(cpm[held+1])-response_baseline
        for method, prediction in values.items():
            expected = np.maximum(expected_raw[method], 0)
            np.testing.assert_allclose(prediction, expected, atol=1e-11, rtol=1e-11)
            diagnostic = diagnostics[method]
            assert diagnostic['renormalized'] is False
            implied = float(np.expm1(prediction).sum())
            np.testing.assert_allclose(implied, diagnostic['impliedCPMSum'], atol=1e-7, rtol=1e-12)
            clipped = int(np.count_nonzero(expected_raw[method] < 0))
            if clipped != diagnostic['clippedFeatures']:
                ambiguous = int(np.count_nonzero(abs(expected_raw[method]) <= 1e-11))
                assert abs(clipped-diagnostic['clippedFeatures']) <= ambiguous
                diagnostic_differences.append(dict(gemgroup=gem, target=target, method=method,
                    nativeClips=diagnostic['clippedFeatures'], referenceClips=clipped, nearZeroFeatures=ambiguous))
            comparisons.append(dict(gemgroup=gem, target=target, method=method,
                maximumPredictionError=float(np.max(abs(prediction-expected))), maximumWeightError=maximum_weight_error))
            for panel, ix in [('allGenes', np.arange(33694)), ('trainingTop1000', top)]:
                scores.append(dict(gemgroup=gem, target=target, method=method, panel=panel,
                    supported=target in annotations, **metrics((prediction-response_baseline)[ix], truth[ix])))
    assert len(comparisons) == 750 and len(scores) == 1500
    summary = []
    for gem in sorted(raw_by_context):
        for panel in ['allGenes', 'trainingTop1000']:
            subset_scores = [x for x in scores if x['gemgroup']==gem and x['panel']==panel and x['supported']]
            means = {method: float(np.mean([x['rmse'] for x in subset_scores if x['method']==method])) for method in METHODS}
            zero = {x['target']:x['rmse'] for x in subset_scores if x['method']=='noChange'}
            mean = {x['target']:x['rmse'] for x in subset_scores if x['method']=='meanSingleResponse'}
            kernel = [x for x in subset_scores if x['method']=='goRidgeFixed']
            summary.append(dict(gemgroup=gem, panel=panel, supportedTargets=len(kernel), meanResponseRMSE=means,
                beatsAllThreeComparators=all(means['goRidgeFixed']<means[m] for m in ['meanSingleResponse','meanSupportedResponse','shuffledGoRidgeFixed']),
                kernelWorseThanNoChange=[x['target'] for x in kernel if x['rmse']>zero[x['target']]],
                kernelWorseThanMean=[x['target'] for x in kernel if x['rmse']>mean[x['target']]]))
    for name, value in [('comparisons', comparisons), ('scores', scores), ('summary', summary),
                        ('diagnostic-differences', diagnostic_differences), ('panels', panels)]:
        write(out/(name+'.json'), value)
    write(out/'checks.json', dict(status='passed-numerical-check-and-complete-scoring', folds=150, comparedVectors=750,
        nativePredictionFreezeSHA256=sha(native/'prediction-freeze.json'), sourceReferenceSHA256=REFERENCE,
        checkerSHA256=sha(Path(__file__)), exactTrainingCounts=True, independentCenteredSklearnKernel=True,
        maximumPredictionError=max(x['maximumPredictionError'] for x in comparisons),
        maximumWeightError=max(x['maximumWeightError'] for x in comparisons),
        clippingClassificationDifferences=len(diagnostic_differences),
        primaryGroupsBeatingAllComparators=sum(x['beatsAllThreeComparators'] for x in summary if x['panel']=='allGenes'),
        packages={n:importlib.metadata.version(n) for n in ['numpy','scikit-learn','scipy']},
        files=[dict(path=p.name, bytes=p.stat().st_size, SHA256=sha(p)) for p in sorted(out.iterdir())]))
    print(json.dumps([x for x in summary if x['panel']=='allGenes'], indent=2))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
