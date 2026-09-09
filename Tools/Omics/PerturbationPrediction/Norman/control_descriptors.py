#!/usr/bin/env python3
"""Control-only sparse gene descriptors and separate held-target prediction/scoring."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from scipy import sparse
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits
from combinations import REFERENCE_SHA256, counts, load, metrics, sha, write
from prepare import SOURCE, column
from unseen_targets import select

METHODS = ('noChange', 'meanSingleResponse', 'controlCorrelationRidge', 'shuffledControlCorrelationRidge')


def identity():
    root = Path(__file__).parent
    return {n: sha(root / n) for n in ['control_descriptors.py', 'CONTROL_DESCRIPTORS_PROTOCOL.md',
                                      'combinations.py', 'unseen_targets.py', 'prepare.py']}


def verify_reference(path):
    if sha(path) != REFERENCE_SHA256:
        raise ValueError('reference fingerprint mismatch')
    return load(path)


def prepare(source_path, reference_path, out):
    ref = verify_reference(reference_path)
    if sha(source_path) != SOURCE['sha256']:
        raise ValueError('source fingerprint mismatch')
    out.mkdir(parents=True, exist_ok=False)
    mask = ref['cell_conditions'] == 'control'
    totals = ref['totals'][mask].astype(float)
    assert np.all(totals > 0) and len(totals) == 11855
    n, p = len(totals), len(ref['feature_ids'])
    cell_map = np.full(len(mask), -1, dtype=np.int64)
    cell_map[mask] = np.arange(n)
    symbols = ref['gene_symbols'].tolist()
    targets = sorted(t for t in ref['conditions'] if t != 'control' and '_' not in t)
    matched = {t: symbols.index(t) for t in targets if symbols.count(t) == 1}
    sums, squares = np.zeros(p), np.zeros(p)
    detected = np.zeros(p, dtype=np.int64)
    raw_totals = np.zeros(n, dtype=np.uint64)
    with h5py.File(source_path) as f:
        np.testing.assert_array_equal(column(f['obs'], '_index'), ref['barcodes'])
        np.testing.assert_array_equal(column(f['obs'], 'perturbation'), ref['cell_conditions'])
        np.testing.assert_array_equal(column(f['var'], 'ensemble_id'), ref['feature_ids'])
        np.testing.assert_array_equal(column(f['var'], '_index'), ref['gene_symbols'])
        x = f['X']
        assert x.attrs['encoding-type'] == 'csc_matrix'
        ptr = x['indptr'][:]
        def read(j):
            start, end = int(ptr[j]), int(ptr[j + 1])
            row = cell_map[x['indices'][start:end]]
            keep = row >= 0
            value = x['data'][start:end][keep].astype(np.float64)
            return row[keep], value
        for j in range(p):
            rows, raw = read(j)
            raw_totals[rows] += raw.astype(np.uint64)
            values = np.log1p(raw / totals[rows] * 1e6)
            sums[j], squares[j], detected[j] = values.sum(), values @ values, len(values)
        np.testing.assert_array_equal(raw_totals, ref['totals'][mask])
        mean = sums / n
        variance = np.maximum(squares / n - mean**2, 0)
        excluded = set(matched.values())
        eligible = [j for j in range(p) if j not in excluded and detected[j] >= 10 and variance[j] > 0]
        landmarks = sorted(eligible, key=lambda j: (-variance[j], ref['feature_ids'][j]))[:2000]
        assert len(landmarks) == 2000
        supported = [t for t in targets if t in matched and detected[matched[t]] >= 10 and variance[matched[t]] > 0]
        target_indices = [matched[t] for t in supported]
        chosen = target_indices + landmarks
        indices, values, offsets = [], [], [0]
        for j in chosen:
            rows, raw = read(j)
            indices.append(rows)
            values.append(np.log1p(raw / totals[rows] * 1e6))
            offsets.append(offsets[-1] + len(rows))
        matrix = sparse.csc_matrix((np.concatenate(values), np.concatenate(indices), np.array(offsets)), shape=(n, len(chosen)))
    a, b = matrix[:, :len(supported)], matrix[:, len(supported):]
    covariance = (a.T @ b).toarray() / n - mean[target_indices, None] * mean[None, landmarks]
    correlation = covariance / np.sqrt(variance[target_indices, None] * variance[None, landmarks])
    assert np.all(np.isfinite(correlation)) and np.max(abs(correlation)) <= 1 + 1e-10
    # Small independent dense oracle: two cell vectors at a time, never a dense
    # control-cell by full-feature or selected-feature matrix.
    errors = []
    for i, j in [(0, 0), (len(supported)//2, 999), (len(supported)-1, 1999)]:
        exact = np.corrcoef(a[:, i].toarray().ravel(), b[:, j].toarray().ravel())[0, 1]
        errors.append(abs(float(exact - correlation[i, j])))
    assert max(errors) < 1e-10
    np.savez_compressed(out / 'descriptors.npz', targetIDs=np.array(supported),
        targetFeatureIDs=ref['feature_ids'][target_indices], landmarkIDs=ref['feature_ids'][landmarks],
        expression=correlation, controlBarcodes=ref['barcodes'][mask], featureIDs=ref['feature_ids'],
        controlMean=mean, controlVariance=variance, controlDetected=detected)
    coverage = [dict(target=t, status='supported' if t in supported else 'unresolved-symbol' if t not in matched else 'insufficient-control-expression',
                     detectedControlCells=int(detected[matched[t]]) if t in matched else None) for t in targets]
    write(out / 'coverage.json', coverage)
    write(out / 'receipt.json', dict(sourceSHA256=SOURCE['sha256'], referenceSHA256=REFERENCE_SHA256,
        implementation=identity(), descriptorSHA256=sha(out / 'descriptors.npz'), coverageSHA256=sha(out / 'coverage.json'),
        controlCells=n, supportedTargets=len(supported), landmarks=2000, sparseSelectedEntries=matrix.nnz,
        controlTotalsExact=True, independentCorrelationMaxAbsoluteError=max(errors), perturbedCellsUsed=0))


def fit_predict(data, descriptors):
    targets = data['trainingTargets']
    raw, ctrl = counts(data['trainingCounts']), counts(data['control'])
    cpm = raw / raw.sum(axis=1, keepdims=True) * 1e6
    baseline = np.log1p(ctrl / ctrl.sum() * 1e6)
    delta = np.log1p(cpm) - baseline
    mean_cpm = (cpm.sum(axis=0) + np.expm1(baseline)) / (len(targets) + 1)
    top = np.array(sorted(range(len(ctrl)), key=lambda i: (-mean_cpm[i], data['featureIDs'][i]))[:1000])
    raw_predictions = dict(noChange=baseline.copy(), meanSingleResponse=baseline + delta.mean(axis=0))
    supported = descriptors['targetIDs'].tolist()
    info = dict(target=data['target'], trainingTargets=targets, excludedConditions=data['excluded'], supported=data['target'] in supported)
    if info['supported']:
        usable = [i for i, t in enumerate(targets) if t in supported]
        x = descriptors['expression'][[supported.index(targets[i]) for i in usable]]
        query = descriptors['expression'][supported.index(data['target'])]
        center, scale = x.mean(axis=0), x.std(axis=0)
        scale[scale == 0] = 1
        x = (x - center) / scale / np.sqrt(x.shape[1])
        query = (query - center) / scale / np.sqrt(x.shape[1])
        weight = np.linalg.solve(x @ x.T + np.eye(len(x)), x @ query)
        errors = []
        for method, y in [(METHODS[2], delta[usable]), (METHODS[3], np.roll(delta[usable], -1, axis=0))]:
            mean = y.mean(axis=0)
            response = mean + weight @ (y - mean)
            independent = Ridge(alpha=1, solver='svd').fit(x, y).predict(query[None])[0]
            errors.append(float(np.max(abs(response - independent))))
            np.testing.assert_allclose(response, independent, rtol=1e-10, atol=1e-10)
            raw_predictions[method] = baseline + response
        info.update(descriptorTrainingTargets=[targets[i] for i in usable], independentMaxAbsoluteError=max(errors))
    predicted = {k: np.maximum(v, 0) for k, v in raw_predictions.items()}
    info['diagnostics'] = {k: dict(clippedFeatures=int(np.count_nonzero(raw_predictions[k] < 0)), impliedCPMSum=float(np.expm1(v).sum())) for k, v in predicted.items()}
    return predicted, baseline, top, info


def predict(reference_path, descriptor_path, out):
    ref = verify_reference(reference_path)
    receipt = json.loads((descriptor_path / 'receipt.json').read_text())
    if sha(descriptor_path / 'descriptors.npz') != receipt['descriptorSHA256'] or receipt['referenceSHA256'] != REFERENCE_SHA256:
        raise ValueError('descriptor provenance mismatch')
    d = load(descriptor_path / 'descriptors.npz')
    np.testing.assert_array_equal(d['featureIDs'], ref['feature_ids'])
    targets = sorted(t for t in ref['conditions'] if t != 'control' and '_' not in t)
    out.mkdir(parents=True, exist_ok=False)
    arrays = {k: [] for k in METHODS}
    infos, baselines, panels, supported = [], [], [], []
    for t in targets:
        data = select(ref, t)
        assert all(t not in x.split('_') for x in data['trainingTargets'])
        predictions, baseline, top, info = fit_predict(data, d)
        for k in METHODS:
            arrays[k].append(predictions.get(k, np.full(len(baseline), np.nan)))
        infos.append(info)
        baselines.append(baseline)
        panels.append(top)
        supported.append(info['supported'])
        print(json.dumps(dict(target=t, supported=info['supported'])), flush=True)
    np.savez_compressed(out / 'predictions.npz', **{k: np.array(v) for k, v in arrays.items()},
                        targetIDs=np.array(targets), featureIDs=ref['feature_ids'], baseline=np.array(baselines),
                        trainingTop1000=np.array(panels), supported=np.array(supported))
    write(out / 'folds.json', infos)
    write(out / 'receipt.json', dict(referenceSHA256=REFERENCE_SHA256, implementation=identity(),
        descriptorReceiptSHA256=sha(descriptor_path / 'receipt.json'), predictionsSHA256=sha(out / 'predictions.npz'),
        foldsSHA256=sha(out / 'folds.json'), independentMaxAbsoluteError=max(i.get('independentMaxAbsoluteError', 0) for i in infos)))


def score(reference_path, prediction_path, out):
    ref = verify_reference(reference_path)
    receipt = json.loads((prediction_path / 'receipt.json').read_text())
    if sha(prediction_path / 'predictions.npz') != receipt['predictionsSHA256']:
        raise ValueError('prediction fingerprint mismatch')
    pred = load(prediction_path / 'predictions.npz')
    np.testing.assert_array_equal(ref['feature_ids'], pred['featureIDs'])
    names = ref['conditions'].tolist()
    results = []
    for i, target in enumerate(pred['targetIDs']):
        raw = counts(ref['counts'][names.index(target)])
        truth = np.log1p(raw / raw.sum() * 1e6) - pred['baseline'][i]
        for method in METHODS:
            if method in METHODS[2:] and not pred['supported'][i]:
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
    write(out / 'receipt.json', dict(referenceSHA256=REFERENCE_SHA256, implementation=identity(),
                                     predictionReceiptSHA256=sha(prediction_path / 'receipt.json')))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['prepare', 'predict', 'score'])
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--source', type=Path)
    p.add_argument('--descriptors', type=Path)
    p.add_argument('--predictions', type=Path)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    with threadpool_limits(limits=1):
        required = {'prepare': a.source, 'predict': a.descriptors, 'score': a.predictions}[a.mode]
        if required is None:
            p.error('mode-specific input is required')
        if a.mode == 'prepare':
            prepare(a.source, a.reference, a.out)
        elif a.mode == 'predict':
            predict(a.reference, a.descriptors, a.out)
        else:
            score(a.reference, a.predictions, a.out)

if __name__ == '__main__':
    main()
