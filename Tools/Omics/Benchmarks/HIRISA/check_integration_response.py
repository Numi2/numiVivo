#!/usr/bin/env python3
"""Stream full-cohort PCs into matched experimental-response diagnostics.

Experimental library assignments are targets; author cell-type predictions are
not ground truth. This is a transductive preservation diagnostic, not prospective
prediction, a biological identity test, or a complete integration benchmark.
"""
import argparse
import hashlib
import itertools
import json
import platform
import resource
import time
from pathlib import Path

import h5py
import ijson
import numpy as np

RECORD = np.dtype([('row', '<u4'), ('column', '<u4'), ('value', '<f8')])


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def blocks(path, n, d, rows):
    assert 1 <= rows <= 8192 and 1 <= d <= 64 and 1 <= n <= 2000000
    assert path.is_file() and not path.is_symlink() and path.stat().st_size == n*d*16
    with path.open('rb') as f:
        for start in range(0, n, rows):
            end = min(n, start+rows)
            raw = f.read((end-start)*d*16)
            assert len(raw) == (end-start)*d*16, 'Truncated score records'
            a = np.frombuffer(raw, dtype=RECORD).reshape(end-start, d)
            assert np.array_equal(a['row'], np.broadcast_to(np.arange(start, end)[:, None], a.shape)), 'Row coordinates'
            assert np.array_equal(a['column'], np.broadcast_to(np.arange(d), a.shape)), 'Column coordinates'
            x = a['value'].copy()
            assert np.isfinite(x).all(), 'Nonfinite score'
            yield start, end, x
        assert not f.read(1)


def one(sample, key):
    assert len(sample[key]) == 1
    return sample[key][0]


def validate_design(design):
    assert design['assignmentUsesExpression'] is False
    samples = design['samples']
    ids = [s['accession'] for s in samples]
    assert len(ids) == len(set(ids)) and 2 <= len(ids) <= 4096
    lookup = dict(zip(ids, samples))
    seen = set()
    for comparison in design['comparisons']:
        key = (comparison['population'], comparison['treatment'])
        assert key not in seen
        seen.add(key)
        donors = [p['donor'] for p in comparison['pairs']]
        assert len(donors) >= 3 and len(donors) == len(set(donors))
        for pair in comparison['pairs']:
            for role, condition in [('control', 'none'), ('treated', comparison['treatment'])]:
                s = lookup[pair[role]]
                assert one(s, 'subject id') == pair['donor']
                assert one(s, 'cell type') == comparison['population']
                assert one(s, 'treatment') == condition
                assert one(s, 'batch id') == pair['batch'] and one(s, 'pool id') == pair['pool']
    return samples, ids


def source_codes(source, metadata, samples, n, rows):
    """Only n uint16 library codes persist; strings/cell dictionaries are tiled."""
    ids = [s['accession'] for s in samples]
    indices = {s: i for i, s in enumerate(ids)}
    codes = np.empty(n, dtype=np.uint16)
    with metadata.open('rb') as f:
        native_samples = list(ijson.items(f, 'samples.item'))
    assert len(native_samples) == len(samples)
    by_id = {s['id']: s for s in native_samples}
    assert set(by_id) == set(ids)
    for s in samples:
        native = by_id[s['accession']]
        assert native['donorID'] == one(s, 'subject id')
        assert native['condition'] == one(s, 'treatment')
        assert native['batchID'] == one(s, 'batch pool')
    with h5py.File(source, 'r') as h, metadata.open('rb') as f:
        assert int(h['X'].attrs['shape'][0]) == n
        fields = {'geo_donor': 'subject id', 'geo_treatment': 'treatment',
                  'geo_enrichment': 'cell type', 'geo_batch_pool': 'batch pool'}
        cells = ijson.items(f, 'cells.item')
        for start in range(0, n, rows):
            end = min(n, start+rows)
            chunk = list(itertools.islice(cells, end-start))
            assert len(chunk) == end-start, 'Metadata row count'
            accessions = h['obs/geo_accession'].asstr()[start:end]
            barcodes = h['obs/_index'].asstr()[start:end]
            assert [c['sampleID'] for c in chunk] == accessions.tolist(), 'Source library order'
            assert [c['barcode'] for c in chunk] == barcodes.tolist(), 'Source barcode order'
            cc = np.array([indices[v] for v in accessions], dtype=np.uint16)
            codes[start:end] = cc
            for field, key in fields.items():
                expected = np.array([one(samples[i], key) for i in cc])
                assert np.array_equal(h['obs/'+field].asstr()[start:end], expected), field
        assert next(cells, None) is None, 'Extra metadata rows'
    assert np.all(np.bincount(codes, minlength=len(samples)) > 0)
    return codes


def moments(path, codes, count, d, rows):
    size = np.bincount(codes, minlength=count)
    sums = np.zeros((count, d))
    products = np.zeros((count, d, d))
    for start, end, x in blocks(path, len(codes), d, rows):
        cc = codes[start:end]
        for k in np.unique(cc):
            y = x[cc == k]
            sums[k] += y.sum(axis=0)
            products[k] += y.T @ y
    means = sums / size[:, None]
    second = products / size[:, None, None]
    assert np.isfinite(means).all() and np.isfinite(second).all()
    return size, means, second


def fit(means, second, selected, targets, ridge):
    """Equal library weights; train-only centering/scaling; unpenalized intercept."""
    assert ridge > 0 and len(selected) == len(targets)
    m = means[selected].mean(axis=0)
    q = second[selected].mean(axis=0)
    variance = np.diag(q) - m*m
    assert variance.min() >= -1e-10 * max(1., float(np.abs(q).max()))
    scale = np.sqrt(np.maximum(variance, 0.))
    scale[scale < 1e-12] = 1.
    gram = (q - np.outer(m, m)) / np.outer(scale, scale)
    target_mean = float(np.mean(targets))
    rhs = ((means[selected] - m) / scale).T @ np.asarray(targets) / len(selected)
    beta = np.linalg.solve(gram + ridge*np.eye(len(m)), rhs)
    weight = beta/scale
    intercept = target_mean - m @ weight
    assert np.isfinite(weight).all() and np.isfinite(intercept)
    return weight, float(intercept)


def evaluate(path, codes, samples, ids, design, d, rows, ridge, library_erasure=False):
    size, means, second = moments(path, codes, len(samples), d, rows)
    original_means = means.copy()
    if library_erasure:
        second = second - np.einsum('si,sj->sij', means, means)
        means = np.zeros_like(means)
    index = {name: i for i, name in enumerate(ids)}
    folds = []
    for c in design['comparisons']:
        for held in c['pairs']:
            train = [p for p in c['pairs'] if p['donor'] != held['donor']]
            selected = [index[p[role]] for p in train for role in ['control', 'treated']]
            targets = [-1., 1.] * len(train)
            w, b = fit(means, second, selected, targets, ridge)
            folds.append(dict(population=c['population'], treatment=c['treatment'], donor=held['donor'],
                              control=index[held['control']], treated=index[held['treated']],
                              weight=w, intercept=b, correct=[0, 0], total=[0, 0],
                              trainingLibraries=[ids[i] for i in selected]))
    by_sample = {}
    for f in folds:
        for role in ['control', 'treated']:
            by_sample.setdefault(f[role], []).append((f, 0 if role == 'control' else 1))
    for start, end, x in blocks(path, len(codes), d, rows):
        cc = codes[start:end]
        if library_erasure:
            x -= original_means[cc]
        for k in np.unique(cc):
            y = x[cc == k]
            for f, target in by_sample.get(int(k), []):
                prediction = (y @ f['weight'] + f['intercept']) >= 0
                f['correct'][target] += int((prediction == target).sum())
                f['total'][target] += len(y)
    for f in folds:
        assert f['total'] == [int(size[f['control']]), int(size[f['treated']])]
        f['recall'] = [a/b for a, b in zip(f['correct'], f['total'])]
        f['balancedAccuracy'] = float(np.mean(f['recall']))
        f['weight'] = f['weight'].tolist()
        f['control'] = ids[f['control']]
        f['treated'] = ids[f['treated']]
    # All libraries contribute to this descriptive donor-associated scatter.
    # Equal library and donor weighting prevents library size from deciding it.
    strata = {}
    for k, s in enumerate(samples):
        key = (one(s, 'cell type'), one(s, 'treatment'))
        strata.setdefault(key, {}).setdefault(one(s, 'subject id'), []).append(k)
    scatter = []
    for (pop, condition), donors in sorted(strata.items()):
        if len(donors) < 2:
            scatter.append(dict(population=pop, condition=condition, status='unavailable-single-donor'))
            continue
        centroids = np.array([means[indices].mean(axis=0) for indices in donors.values()])
        center = centroids.mean(axis=0)
        between = float(np.mean(np.sum((centroids-center)**2, axis=1)))
        total = float(np.mean([np.trace(second[indices].mean(axis=0)) for indices in donors.values()]) - center @ center)
        scatter.append(dict(population=pop, condition=condition, donors=len(donors),
                            betweenDonorSquaredDistance=between, totalVariance=max(0., total),
                            fraction=between/total if total > 1e-12 else None,
                            status='measured' if total > 1e-12 else 'unavailable-zero-variance'))
    return dict(cells=len(codes), libraryCells={key: int(n) for key, n in zip(ids, size)},
                libraryMeans=means.tolist(), folds=folds, conditionalDonorScatter=scatter,
                classificationCells=int(sum(size[list(by_sample)])),
                cellsOutsideMatchedClassifiers=int(len(codes)-sum(size[list(by_sample)])))


def compare(baseline, candidate, design, ids, margins):
    assert baseline['libraryCells'] == candidate['libraryCells']
    x = np.asarray(baseline['libraryMeans']); y = np.asarray(candidate['libraryMeans'])
    index = {name: i for i, name in enumerate(ids)}
    results = []
    for c in design['comparisons']:
        key = (c['population'], c['treatment'])
        base_folds = [f for f in baseline['folds'] if (f['population'], f['treatment']) == key]
        new_folds = [f for f in candidate['folds'] if (f['population'], f['treatment']) == key]
        assert [f['donor'] for f in base_folds] == [f['donor'] for f in new_folds]
        loss = [a['balancedAccuracy']-b['balancedAccuracy'] for a, b in zip(base_folds, new_folds)]
        responses = []
        for pair in c['pairs']:
            a = x[index[pair['treated']]]-x[index[pair['control']]]
            b = y[index[pair['treated']]]-y[index[pair['control']]]
            norm = float(np.linalg.norm(a)); other = float(np.linalg.norm(b))
            if norm <= margins['minimumResponseNorm']:
                responses.append(dict(donor=pair['donor'], status='unavailable-negligible-baseline-response'))
            else:
                responses.append(dict(donor=pair['donor'], status='measured', baselineNorm=norm,
                                      candidateNorm=other, cosine=float(a@b/(norm*other)) if other > 1e-12 else None,
                                      relativeDrift=float(np.linalg.norm(b-a)/norm),
                                      projectedGain=float(a@b/(norm*norm))))
        complete = all(v['status'] == 'measured' for v in responses)
        mean_drift = float(np.mean([v['relativeDrift'] for v in responses])) if complete else None
        results.append(dict(population=key[0], treatment=key[1], folds=len(loss),
                            meanBalancedAccuracyLoss=float(np.mean(loss)), maximumFoldLoss=max(loss),
                            baselineBalancedAccuracy=float(np.mean([f['balancedAccuracy'] for f in base_folds])),
                            candidateBalancedAccuracy=float(np.mean([f['balancedAccuracy'] for f in new_folds])),
                            pairedResponses=responses, meanRelativeResponseDrift=mean_drift,
                            gates=dict(meanClassificationPreserved=float(np.mean(loss)) <= margins['maximumMeanAccuracyLoss'],
                                       everyFoldPreserved=max(loss) <= margins['maximumFoldAccuracyLoss'],
                                       completeResponsePairs=complete,
                                       meanResponsePreserved=mean_drift is not None and mean_drift <= margins['maximumMeanResponseDrift'])))
    return results


def validate_baseline(baseline, bindings, evaluator_sha, protocol):
    assert baseline['status'] == 'measured' and not baseline['libraryMeanErasure']
    assert 'comparisons' not in baseline, 'A corrected result cannot be a baseline'
    for name in ['source', 'metadata', 'design', 'protocol']:
        assert baseline['bindings'][name] == bindings[name], ('Baseline identity', name)
    assert baseline['bindings']['scores']['SHA256'] == protocol['baselineScoresSHA256']
    assert baseline['evaluatorSHA256'] == evaluator_sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'metadata', 'design', 'protocol', 'scores', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--baseline', type=Path)
    p.add_argument('--library-erasure', action='store_true')
    a = p.parse_args(); assert not a.out.exists()
    began = time.time(); protocol = json.loads(a.protocol.read_text())
    inputs = {name: getattr(a, name) for name in ['source', 'metadata', 'design', 'protocol', 'scores']}
    if a.baseline:
        inputs['baseline'] = a.baseline
    bindings = {name: dict(bytes=path.stat().st_size, SHA256=sha(path)) for name, path in inputs.items()}
    for name in ['source', 'metadata', 'design']:
        assert bindings[name]['SHA256'] == protocol[name+'SHA256'], ('Frozen input', name)
    assert 0 < protocol['ridge'] and 1 <= protocol['rowBatch'] <= 8192
    design = json.loads(a.design.read_text()); samples, ids = validate_design(design)
    assert len(samples) == protocol['libraries']
    codes = source_codes(a.source, a.metadata, samples, protocol['cells'], protocol['rowBatch'])
    print('All source identities and experimental library assignments checked', flush=True)
    result = evaluate(a.scores, codes, samples, ids, design, protocol['components'], protocol['rowBatch'], protocol['ridge'], a.library_erasure)
    result.update(status='measured', libraryMeanErasure=a.library_erasure, bindings=bindings,
                  protocol=protocol, evaluatorSHA256=sha(Path(__file__)),
                  versions=dict(numpy=np.__version__, h5py=h5py.__version__, ijson=ijson.__version__),
                  scope=protocol['scope'])
    if a.baseline:
        baseline = json.loads(a.baseline.read_text())
        validate_baseline(baseline, bindings, result['evaluatorSHA256'], protocol)
        result['comparisons'] = compare(baseline, result, design, ids, protocol['margins'])
        result['allResponseDiagnosticGatesPassed'] = all(all(c['gates'].values()) for c in result['comparisons'])
    else:
        assert not a.library_erasure and bindings['scores']['SHA256'] == protocol['baselineScoresSHA256']
    assert all(path.stat().st_size == bindings[name]['bytes'] and sha(path) == bindings[name]['SHA256'] for name, path in inputs.items()), 'Input changed'
    result['seconds'] = time.time()-began
    result['platform'] = platform.platform()
    result['maximumResidentBytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if platform.system() == 'Darwin' else 1024)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('x') as f:
        json.dump(result, f, indent=2, sort_keys=True, allow_nan=False); f.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'cells', 'classificationCells', 'cellsOutsideMatchedClassifiers', 'seconds', 'maximumResidentBytes']}))


if __name__ == '__main__':
    main()
