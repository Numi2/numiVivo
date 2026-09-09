#!/usr/bin/env python3
"""LOTO co-response reference; prediction and scoring are separate commands."""
import argparse
import json
from pathlib import Path
import numpy as np
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits
from combinations import REFERENCE_SHA256, counts, load, metrics, sha, write

METHODS = ('noChange', 'meanSingleResponse', 'coResponseRidge', 'shuffledCoResponseRidge')


def select(ref, target):
    names = ref['conditions'].tolist()
    singles = sorted(n for n in names if n != 'control' and '_' not in n)
    if target not in singles or len(set(names)) != len(names):
        raise ValueError('unique conditions and known single target required')
    train = [n for n in singles if n != target]
    excluded = sorted(n for n in names if target in n.split('_'))
    return dict(target=target, trainingTargets=train, excluded=excluded,
                featureIDs=ref['feature_ids'], symbols=ref['gene_symbols'],
                control=ref['counts'][names.index('control')].copy(),
                trainingCounts=ref['counts'][[names.index(n) for n in train]].copy())


def fit_predict(data):
    raw, ctrl = counts(data['trainingCounts']), counts(data['control'])
    features, symbols, targets = data['featureIDs'], data['symbols'].tolist(), data['trainingTargets']
    if raw.shape != (len(targets), len(features)) or ctrl.shape != (len(features),):
        raise ValueError('matrix axes disagree')
    if len(set(features)) != len(features) or len(symbols) != len(features):
        raise ValueError('invalid feature identities')
    cpm = raw / raw.sum(axis=1, keepdims=True) * 1e6
    baseline = np.log1p(ctrl / ctrl.sum() * 1e6)
    delta = np.log1p(cpm) - baseline
    mean_cpm = (cpm.sum(axis=0) + np.expm1(baseline)) / (len(targets) + 1)
    top = np.array(sorted(range(len(features)), key=lambda i: (-mean_cpm[i], features[i]))[:1000])
    unbounded = dict(noChange=baseline.copy(), meanSingleResponse=baseline + delta.mean(axis=0))
    info = dict(target=data['target'], trainingTargets=targets, excludedConditions=data['excluded'],
                descriptorStatus='unsupported-symbol', independentMaxAbsoluteError=None)
    if symbols.count(data['target']) == 1:
        usable = [i for i, t in enumerate(targets) if symbols.count(t) == 1]
        x = delta[:, [symbols.index(targets[i]) for i in usable]].T
        query = delta[:, symbols.index(data['target'])]
        center, scale = x.mean(axis=0), x.std(axis=0)
        scale[scale == 0] = 1
        x = (x - center) / scale / np.sqrt(x.shape[1])
        query = (query - center) / scale / np.sqrt(x.shape[1])
        y = delta[usable]
        weight = np.linalg.solve(x @ x.T + np.eye(len(x)), x @ query)
        errors = []
        for method, output in [('coResponseRidge', y), ('shuffledCoResponseRidge', np.roll(y, -1, axis=0))]:
            mean = output.mean(axis=0)
            response = mean + weight @ (output - mean)
            independent = Ridge(alpha=1, solver='svd').fit(x, output).predict(query[None])[0]
            errors.append(float(np.max(abs(response - independent))))
            np.testing.assert_allclose(response, independent, rtol=1e-10, atol=1e-10)
            unbounded[method] = baseline + response
        info.update(descriptorStatus='supported', descriptorTrainingTargets=[targets[i] for i in usable],
                    descriptorCoordinates=targets, independentMaxAbsoluteError=max(errors))
    predictions = {k: np.maximum(v, 0) for k, v in unbounded.items()}
    info['diagnostics'] = {k: dict(clippedFeatures=int(np.count_nonzero(unbounded[k] < 0)),
                                  impliedCPMSum=float(np.expm1(v).sum()), renormalized=False)
                           for k, v in predictions.items()}
    return predictions, baseline, top, info


def predict(reference_path, out):
    if sha(reference_path) != REFERENCE_SHA256:
        raise ValueError('reference hash differs from full qualified source')
    ref = load(reference_path)
    names = ref['conditions'].tolist()
    targets = sorted(n for n in names if n != 'control' and '_' not in n)
    if len(targets) != 105:
        raise ValueError('expected 105 folds')
    out.mkdir(parents=True, exist_ok=False)
    records, baseline_rows, panels, supported = [], [], [], []
    arrays = {method: [] for method in METHODS}
    for target in targets:
        data = select(ref, target)
        changed = dict(ref, counts=ref['counts'].copy())
        changed['counts'][[names.index(n) for n in data['excluded']]] = 1
        mutated = select(changed, target)
        for key in ('control', 'trainingCounts', 'featureIDs', 'symbols'):
            np.testing.assert_array_equal(data[key], mutated[key])
        assert all(target not in n.split('_') for n in data['trainingTargets'])
        predictions, baseline, top, info = fit_predict(data)
        info['excludedOutcomeMutationInputsExact'] = True
        records.append(info)
        baseline_rows.append(baseline)
        panels.append(top)
        supported.append('coResponseRidge' in predictions)
        for method in METHODS:
            arrays[method].append(predictions.get(method, np.full(len(baseline), np.nan)))
        print(json.dumps(dict(target=target, status=info['descriptorStatus'])), flush=True)
    np.savez_compressed(out / 'predictions.npz', **{k: np.array(v) for k, v in arrays.items()},
                        targetIDs=np.array(targets), featureIDs=ref['feature_ids'], baseline=np.array(baseline_rows),
                        trainingTop1000=np.array(panels), supported=np.array(supported))
    write(out / 'folds.json', records)
    write(out / 'receipt.json', dict(referenceSHA256=REFERENCE_SHA256, folds=len(targets),
        protocolSHA256=sha(Path(__file__).with_name('UNSEEN_TARGETS_PROTOCOL.md')),
        implementationSHA256=sha(Path(__file__)), predictionsSHA256=sha(out / 'predictions.npz'),
        foldsSHA256=sha(out / 'folds.json'), outcomesReadByFitter=False,
        independentMaxAbsoluteError=max(r['independentMaxAbsoluteError'] or 0 for r in records)))


def score(reference_path, prediction_path, out):
    receipt = json.loads((prediction_path / 'receipt.json').read_text())
    if sha(reference_path) != receipt['referenceSHA256'] or sha(prediction_path / 'predictions.npz') != receipt['predictionsSHA256']:
        raise ValueError('input fingerprint mismatch')
    ref, pred = load(reference_path), load(prediction_path / 'predictions.npz')
    np.testing.assert_array_equal(ref['feature_ids'], pred['featureIDs'])
    names = ref['conditions'].tolist()
    results = []
    for i, target in enumerate(pred['targetIDs']):
        raw = counts(ref['counts'][names.index(target)])
        truth = np.log1p(raw / raw.sum() * 1e6) - pred['baseline'][i]
        for method in METHODS:
            if not pred['supported'][i] and method in METHODS[2:]:
                continue
            for panel, ix in [('allGenes', np.arange(len(truth))), ('trainingTop1000', pred['trainingTop1000'][i])]:
                result = metrics((pred[method][i] - pred['baseline'][i])[ix], truth[ix])
                results.append(dict(target=str(target), supported=bool(pred['supported'][i]), method=method, panel=panel, **result))
    summary = []
    for cohort in ('allAvailable', 'matchedSupported'):
        for method in METHODS:
            for panel in ('allGenes', 'trainingTop1000'):
                rows = [r for r in results if r['method'] == method and r['panel'] == panel and (cohort == 'allAvailable' or r['supported'])]
                zero = {r['target']: r['rmse'] for r in results if r['method'] == 'noChange' and r['panel'] == panel}
                summary.append(dict(cohort=cohort, method=method, panel=panel, targets=len(rows),
                    meanResponseRMSE=float(np.mean([r['rmse'] for r in rows])),
                    worseThanNoChange=[r['target'] for r in rows if r['rmse'] > zero[r['target']]]))
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'results.json', results)
    write(out / 'summary.json', summary)
    write(out / 'receipt.json', dict(predictionReceiptSHA256=sha(prediction_path / 'receipt.json'),
                                     referenceSHA256=sha(reference_path), implementationSHA256=sha(Path(__file__))))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['predict', 'score'])
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--predictions', type=Path)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    with threadpool_limits(limits=1):
        if a.mode == 'predict':
            predict(a.reference, a.out)
        else:
            if a.predictions is None:
                p.error('score requires --predictions')
            score(a.reference, a.predictions, a.out)

if __name__ == '__main__':
    main()
