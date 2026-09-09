#!/usr/bin/env python3
"""Separate training, prediction and scoring stages for the full Norman pair holdout."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression

REFERENCE_SHA256 = 'ec22f61196cf51f7c1a62c681e4076bd9ea269429f1423dfd65e21c11a7ce8bf'
METHODS = ('noChange', 'meanSingleResponse', 'additiveLogResponse',
           'meanConstituentLogResponse', 'additiveCPMResponse', 'shuffledAdditiveLogResponse')


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def counts(value):
    original = np.asarray(value)
    if original.dtype.kind not in 'uif' or np.any(original > 2**53):
        raise ValueError('count type or integer range is unsupported')
    value = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(value) & (value >= 0) & (value == np.floor(value)) & (value <= 2**53)):
        raise ValueError('counts must be finite nonnegative exact integers')
    if not np.all(value.sum(axis=-1) > 0):
        raise ValueError('empty count library')
    return value


def select_inputs(reference):
    names = reference['conditions'].tolist()
    singles = sorted(n for n in names if n != 'control' and '_' not in n)
    pairs = sorted(n for n in names if '_' in n)
    assert len(names) == len(set(names)) == 237 and len(singles) == 105 and len(pairs) == 131
    queries = [dict(id=n, targets=n.split('_')) for n in pairs]
    assert all(len(q['targets']) == 2 and set(q['targets']) <= set(singles) for q in queries)
    training = dict(featureIDs=reference['feature_ids'], targetIDs=np.array(singles),
                    controlCounts=reference['counts'][names.index('control')],
                    singleCounts=reference['counts'][[names.index(n) for n in singles]])
    truth = dict(featureIDs=reference['feature_ids'], queryIDs=np.array(pairs),
                 counts=reference['counts'][[names.index(n) for n in pairs]])
    return training, queries, truth


def split(reference_path, destination):
    assert sha(reference_path) == REFERENCE_SHA256, 'reference archive differs from full native-qualified input'
    reference = load(reference_path)
    training, queries, truth = select_inputs(reference)
    destination.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(destination / 'training.npz', **training)
    np.savez_compressed(destination / 'sealed-outcomes.npz', **truth)
    write(destination / 'queries.json', queries)
    write(destination / 'inputs.json', dict(referenceSHA256=REFERENCE_SHA256,
                                           trainingSHA256=sha(destination / 'training.npz'),
                                           querySHA256=sha(destination / 'queries.json'),
                                           sealedOutcomeSHA256=sha(destination / 'sealed-outcomes.npz')))


def fit(training):
    control, singles = counts(training['controlCounts']), counts(training['singleCounts'])
    features, targets = training['featureIDs'], training['targetIDs']
    if control.ndim != 1 or singles.shape != (len(targets), len(features)) or len(control) != len(features):
        raise ValueError('training dimensions disagree')
    if len(targets) < 2 or len(set(targets)) != len(targets) or len(set(features)) != len(features):
        raise ValueError('unique training targets and features required')
    order = np.argsort(targets)
    targets, singles = targets[order], singles[order]
    control_cpm = control / control.sum() * 1e6
    single_cpm = singles / singles.sum(axis=1, keepdims=True) * 1e6
    delta = np.log1p(single_cpm) - np.log1p(control_cpm)
    all_training = np.vstack([control, singles])
    selected = np.flatnonzero((all_training.sum(axis=0) >= 10) & ((all_training > 0).sum(axis=0) >= 2))
    mean_cpm = np.vstack([control_cpm, single_cpm]).mean(axis=0)
    top = np.array(sorted(range(len(features)), key=lambda i: (-mean_cpm[i], features[i]))[:1000])
    # Independent least-squares implementation of the single-target identity
    # design; paired outcomes are not arguments to this function.
    linear = LinearRegression(fit_intercept=False).fit(np.eye(len(targets)), delta)
    np.testing.assert_allclose(linear.coef_.T, delta, rtol=1e-12, atol=1e-12)
    return dict(featureIDs=features, targetIDs=targets, controlCPM=control_cpm,
                singleLogResponse=delta, singleCPMResponse=single_cpm - control_cpm,
                meanSingleResponse=delta.mean(axis=0), trainingExpressed=selected,
                trainingTop1000Expression=top), linear


def predict_one(model, targets):
    known = model['targetIDs'].tolist()
    if len(targets) != 2 or len(set(targets)) != 2 or any(t not in known for t in targets):
        raise ValueError('query requires two distinct observed training targets')
    a, b = (known.index(t) for t in targets)
    baseline = np.log1p(model['controlCPM'])
    additive = model['singleLogResponse'][a] + model['singleLogResponse'][b]
    cpm_raw = model['controlCPM'] + (model['singleCPMResponse'][a] + model['singleCPMResponse'][b])
    raw = dict(noChange=baseline, meanSingleResponse=baseline + model['meanSingleResponse'],
               additiveLogResponse=baseline + additive,
               meanConstituentLogResponse=baseline + additive / 2,
               additiveCPMResponse=np.log1p(np.maximum(cpm_raw, 0)),
               shuffledAdditiveLogResponse=baseline + (model['singleLogResponse'][(a + 1) % len(known)]
                                                      + model['singleLogResponse'][(b + 1) % len(known)]))
    predicted = {name: np.maximum(value, 0) for name, value in raw.items()}
    top = np.array(sorted(range(len(baseline)), key=lambda i: (-abs(additive[i]), model['featureIDs'][i]))[:200])
    diagnostics = {name: dict(impliedCPMSum=float(np.expm1(predicted[name]).sum()),
                              clippedFeatures=int(np.count_nonzero(cpm_raw < 0) if name == 'additiveCPMResponse'
                                                  else np.count_nonzero(raw[name] < 0)),
                              renormalized=False) for name in METHODS}
    return predicted, top, diagnostics


def fit_predict(training_path, query_path, destination):
    training = load(training_path)
    queries = json.loads(query_path.read_text())
    if not queries or len({q['id'] for q in queries}) != len(queries):
        raise ValueError('nonempty unique query identities required')
    model, independent = fit(training)
    destination.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(destination / 'model.npz', **model)
    frozen = load(destination / 'model.npz')
    predictions = {method: [] for method in METHODS}
    panels, diagnostics = [], []
    max_linear_error = 0.0
    for query in queries:
        values, top, info = predict_one(model, query['targets'])
        replay, replay_top, replay_info = predict_one(frozen, query['targets'])
        reverse, reverse_top, reverse_info = predict_one(frozen, list(reversed(query['targets'])))
        for method in METHODS:
            np.testing.assert_array_equal(values[method], replay[method])
            np.testing.assert_array_equal(values[method], reverse[method])
            predictions[method].append(values[method])
        np.testing.assert_array_equal(top, replay_top)
        np.testing.assert_array_equal(top, reverse_top)
        assert info == replay_info == reverse_info
        design = np.array([float(t in query['targets']) for t in model['targetIDs']])[None, :]
        expected = np.maximum(np.log1p(model['controlCPM']) + independent.predict(design)[0], 0)
        error = float(np.max(abs(expected - values['additiveLogResponse'])))
        max_linear_error = max(max_linear_error, error)
        np.testing.assert_allclose(expected, values['additiveLogResponse'], rtol=1e-12, atol=1e-12)
        panels.append(top)
        diagnostics.append(dict(query=query['id'], models=info))
    for method, values in predictions.items():
        # Two baselines intentionally share the same prediction for every pair.
        shared = method in ('noChange', 'meanSingleResponse')
        if shared:
            assert all(np.array_equal(values[0], v) for v in values)
        np.savez_compressed(destination / (method + '.npz'),
                            expression=np.array(values[:1] if shared else values),
                            queryRows=np.zeros(len(queries), dtype=np.int64) if shared else np.arange(len(queries)))
    np.savez_compressed(destination / 'panels.npz', queryTop200Response=np.array(panels))
    write(destination / 'queries.json', queries)
    write(destination / 'diagnostics.json', diagnostics)
    protocol = Path(__file__).with_name('COMBINATIONS_PROTOCOL.md')
    write(destination / 'receipt.json', dict(schemaVersion=1, trainingSHA256=sha(training_path),
                                            querySHA256=sha(query_path), protocolSHA256=sha(protocol),
                                            implementationSHA256=sha(Path(__file__)),
                                            independentLinearMaxAbsoluteError=max_linear_error,
                                            frozenReloadExact=True, targetOrderExact=True,
                                            packages={n: importlib.metadata.version(n) for n in ['numpy', 'scikit-learn']},
                                            files={p.name: sha(p) for p in sorted(destination.iterdir()) if p.is_file()}))


def metrics(predicted, observed):
    error = predicted - observed
    active = observed != 0
    denom = float(np.sum(observed**2))
    return dict(features=len(observed), rmse=float(np.sqrt(np.mean(error**2))),
                mae=float(np.mean(abs(error))),
                pearson=float(np.corrcoef(predicted, observed)[0, 1]) if np.std(predicted) > 0 and np.std(observed) > 0 else None,
                nonzeroTruthFeatures=int(active.sum()),
                signAgreement=float(np.mean(np.sign(predicted[active]) == np.sign(observed[active]))) if active.any() else None,
                explainedSumSquaresVsZero=float(1 - np.sum(error**2) / denom) if denom > 0 else None)


def score(truth_path, prediction_path, destination):
    receipt = json.loads((prediction_path / 'receipt.json').read_text())
    assert all(sha(prediction_path / name) == digest for name, digest in receipt['files'].items())
    truth, model = load(truth_path), load(prediction_path / 'model.npz')
    raw = counts(truth['counts'])
    queries = json.loads((prediction_path / 'queries.json').read_text())
    assert truth['queryIDs'].tolist() == [q['id'] for q in queries]
    np.testing.assert_array_equal(truth['featureIDs'], model['featureIDs'])
    observed = np.log1p(raw / raw.sum(axis=1, keepdims=True) * 1e6)
    baseline = np.log1p(model['controlCPM'])
    saved = {name: load(prediction_path / (name + '.npz')) for name in METHODS}
    panel_file = load(prediction_path / 'panels.npz')
    results = []
    for i, query in enumerate(queries):
        panels = dict(allGenes=np.arange(len(baseline)), trainingExpressed=model['trainingExpressed'],
                      trainingTop1000Expression=model['trainingTop1000Expression'],
                      queryTop200Response=panel_file['queryTop200Response'][i])
        result = dict(query=query['id'], targets=query['targets'], models={})
        for method in METHODS:
            predicted = saved[method]['expression'][saved[method]['queryRows'][i]]
            result['models'][method] = {
                panel: dict(response=metrics((predicted - baseline)[ix], (observed[i] - baseline)[ix]),
                            expression=metrics(predicted[ix], observed[i, ix])) for panel, ix in panels.items()}
        results.append(result)
    summary = {}
    for method in METHODS:
        summary[method] = {}
        for panel in results[0]['models'][method]:
            values = [r['models'][method][panel]['response'] for r in results]
            summary[method][panel] = dict(meanResponseRMSE=float(np.mean([v['rmse'] for v in values])),
                                          medianResponseRMSE=float(np.median([v['rmse'] for v in values])),
                                          meanResponsePearson=float(np.mean([v['pearson'] for v in values if v['pearson'] is not None]))
                                          if any(v['pearson'] is not None for v in values) else None,
                                          worseThanNoChange=[r['query'] for r in results if r['models'][method][panel]['response']['rmse']
                                                             > r['models']['noChange'][panel]['response']['rmse']],
                                          betterThanShuffled=sum(r['models'][method][panel]['response']['rmse']
                                                                 < r['models']['shuffledAdditiveLogResponse'][panel]['response']['rmse'] for r in results))
    destination.mkdir(parents=True, exist_ok=False)
    write(destination / 'results.json', results)
    write(destination / 'summary.json', summary)
    write(destination / 'receipt.json', dict(truthSHA256=sha(truth_path), predictionReceiptSHA256=sha(prediction_path / 'receipt.json'),
                                            queries=len(queries), features=len(baseline), allPairsScored=True,
                                            statement='External pooled-condition baselines only. No native predictor, unseen-target, interaction, donor, Bayesian or mechanistic qualification.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    s = commands.add_parser('split')
    s.add_argument('--reference', type=Path, required=True)
    p = commands.add_parser('predict')
    p.add_argument('--training', type=Path, required=True)
    p.add_argument('--queries', type=Path, required=True)
    e = commands.add_parser('score')
    e.add_argument('--truth', type=Path, required=True)
    e.add_argument('--predictions', type=Path, required=True)
    for command in (s, p, e):
        command.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'split':
        split(args.reference, args.out)
    elif args.command == 'predict':
        fit_predict(args.training, args.queries, args.out)
    else:
        score(args.truth, args.predictions, args.out)
    print(json.dumps(dict(status='passed', stage=args.command)))


if __name__ == '__main__':
    main()
