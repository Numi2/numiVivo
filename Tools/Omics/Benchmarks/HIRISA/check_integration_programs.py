#!/usr/bin/env python3
"""Measure held-out donor program gradients represented by complete-cohort PCs."""
import argparse
import json
from pathlib import Path
import platform
import resource
import time
import numpy as np
from check_integration_response import blocks, sha, source_codes, validate_design


def moments(path, programs, codes, count, d, rows):
    p = programs.shape[1]
    number = np.zeros(count, dtype=np.int64)
    sx = np.zeros((count, d)); sxx = np.zeros((count, d, d))
    sy = np.zeros((count, p)); syy = np.zeros((count, p)); sxy = np.zeros((count, d, p))
    for start, end, x in blocks(path, len(codes), d, rows):
        y = np.asarray(programs[start:end]); cc = codes[start:end]
        assert not np.isinf(y).any()
        valid = np.isfinite(y).all(axis=1)
        assert np.array_equal(np.isnan(y).all(axis=1), ~valid), 'Partially missing program row'
        for k in np.unique(cc):
            selected = (cc == k) & valid
            a = x[selected]; b = y[selected]
            number[k] += len(a); sx[k] += a.sum(axis=0); sxx[k] += a.T @ a
            sy[k] += b.sum(axis=0); syy[k] += (b*b).sum(axis=0); sxy[k] += a.T @ b
    denominator = np.maximum(number, 1)
    return dict(number=number, x=sx/denominator[:, None], xx=sxx/denominator[:, None, None],
                y=sy/denominator[:, None], yy=syy/denominator[:, None], xy=sxy/denominator[:, None, None])


def fit(stats, selected, ridge):
    assert ridge > 0 and np.all(stats['number'][selected] > 0)
    x = stats['x'][selected].mean(axis=0); xx = stats['xx'][selected].mean(axis=0)
    y = stats['y'][selected].mean(axis=0); xy = stats['xy'][selected].mean(axis=0)
    variance = np.diag(xx)-x*x
    assert variance.min() >= -1e-10*max(1., float(np.abs(xx).max()))
    scale = np.sqrt(np.maximum(variance, 0)); scale[scale < 1e-12] = 1
    gram = (xx-np.outer(x, x))/np.outer(scale, scale)
    rhs = (xy-np.outer(x, y))/scale[:, None]
    weight = np.linalg.solve(gram+ridge*np.eye(len(x)), rhs)/scale[:, None]
    intercept = y-x@weight
    assert np.isfinite(weight).all() and np.isfinite(intercept).all()
    return weight, intercept, y


def library_metrics(stats, index, weight, intercept, training_mean, minimum_variance):
    x = stats['x'][index]; y = stats['y'][index]
    covxx = stats['xx'][index]-np.outer(x, x)
    covxy = stats['xy'][index]-np.outer(x, y)
    truth_var = stats['yy'][index]-y*y
    prediction_var = np.sum(weight*(covxx@weight), axis=0)
    covariance = np.sum(weight*covxy, axis=0)
    assert truth_var.min() >= -1e-10 and prediction_var.min() >= -1e-10
    truth_var = np.maximum(truth_var, 0); prediction_var = np.maximum(prediction_var, 0)
    error = truth_var+prediction_var-2*covariance+(x@weight+intercept-y)**2
    assert error.min() >= -1e-10
    denominator = truth_var+(training_mean-y)**2
    result = []
    for j in range(len(y)):
        available = stats['number'][index] > 0 and truth_var[j] > minimum_variance
        result.append(dict(status='measured' if available else 'unavailable-program-variance',
                           cells=int(stats['number'][index]), meanProgram=float(y[j]),
                           programVariance=float(truth_var[j]), predictedVariance=float(prediction_var[j]),
                           mse=float(max(error[j], 0)),
                           withinLibraryR2=float((2*covariance[j]-prediction_var[j])/truth_var[j]) if available else None,
                           programPearson=float(covariance[j]/np.sqrt(truth_var[j]*prediction_var[j])) if available and prediction_var[j] > minimum_variance else None,
                           r2VsTrainingMean=float(1-max(error[j], 0)/denominator[j]) if denominator[j] > minimum_variance else None))
    return result


def evaluate(stats, codes, ids, design, protocol, library_mean_only=False):
    if library_mean_only:
        stats = {k: v.copy() for k, v in stats.items()}
        stats['xx'] = np.einsum('si,sj->sij', stats['x'], stats['x'])
        stats['xy'] = np.einsum('si,sj->sij', stats['x'], stats['y'])
    index = {name: i for i, name in enumerate(ids)}
    folds = []
    for comparison in design['comparisons']:
        for held in comparison['pairs']:
            training = [p for p in comparison['pairs'] if p['donor'] != held['donor']]
            selected = [index[p[role]] for p in training for role in ['control', 'treated']]
            query = [index[held[role]] for role in ['control', 'treated']]
            fold = dict(population=comparison['population'], treatment=comparison['treatment'], donor=held['donor'],
                        trainingLibraries=[ids[i] for i in selected], evaluationLibraries=[ids[i] for i in query])
            if np.any(stats['number'][selected+query] == 0):
                fold.update(status='unavailable-empty-library', measurements=[dict(programID=p, withinLibraryR2=None) for p in protocol['programs']])
            else:
                weight, intercept, mean = fit(stats, selected, protocol['ridge'])
                metrics = [library_metrics(stats, k, weight, intercept, mean, protocol['margins']['minimumProgramVariance']) for k in query]
                measurements = []
                for j, program in enumerate(protocol['programs']):
                    values = [m[j]['withinLibraryR2'] for m in metrics]
                    measurements.append(dict(programID=program, libraries=[m[j] for m in metrics],
                                             withinLibraryR2=float(np.mean(values)) if all(v is not None for v in values) else None))
                fold.update(status='measured', weight=weight.tolist(), intercept=intercept.tolist(), measurements=measurements)
            folds.append(fold)
    matched = sorted({index[name] for fold in folds for name in fold['evaluationLibraries']})
    original_number = np.bincount(codes, minlength=len(ids))
    return dict(cells=len(codes), matchedOriginalCells=int(original_number[matched].sum()),
                matchedProgramTargetCells=int(stats['number'][matched].sum()),
                cellsOutsideMatchedFolds=int(len(codes)-original_number[matched].sum()),
                originalLibraryCells=original_number.tolist(),
                programTargetLibraryCells=stats['number'].tolist(), unavailableEmptyCellTargets=int(len(codes)-stats['number'].sum()),
                moments={k: v.tolist() for k, v in stats.items()}, folds=folds)


def compare(baseline, result, protocol):
    assert baseline['originalLibraryCells'] == result['originalLibraryCells']
    assert baseline['programTargetLibraryCells'] == result['programTargetLibraryCells']
    assert len(baseline['folds']) == len(result['folds'])
    groups = {}
    for first, second in zip(baseline['folds'], result['folds']):
        for key in ['population', 'treatment', 'donor', 'trainingLibraries', 'evaluationLibraries']:
            assert first[key] == second[key]
        assert len(first['measurements']) == len(second['measurements']) == len(protocol['programs'])
        for a, b in zip(first['measurements'], second['measurements']):
            assert a['programID'] == b['programID']
            groups.setdefault((first['population'], first['treatment'], a['programID']), []).append(
                dict(donor=first['donor'], baseline=a['withinLibraryR2'], candidate=b['withinLibraryR2']))
    result = []
    for (population, treatment, program), folds in groups.items():
        complete = all(f['baseline'] is not None and f['candidate'] is not None for f in folds)
        baseline_mean = float(np.mean([f['baseline'] for f in folds])) if complete else None
        losses = [f['baseline']-f['candidate'] for f in folds] if complete else []
        sensitive = complete and baseline_mean > protocol['margins']['minimumErasureDetection']
        result.append(dict(population=population, treatment=treatment, programID=program, folds=folds,
                           status='measured' if complete else 'unavailable-folds',
                           baselineMeanWithinLibraryR2=baseline_mean,
                           meanWithinLibraryR2Loss=float(np.mean(losses)) if complete else None,
                           maximumFoldWithinLibraryR2Loss=max(losses) if complete else None,
                           controlSensitive=sensitive,
                           gates=dict(completeFolds=complete, erasureControlSensitive=sensitive,
                                      meanPreserved=complete and float(np.mean(losses)) <= protocol['margins']['maximumMeanWithinLibraryR2Loss'],
                                      everyFoldPreserved=complete and max(losses) <= protocol['margins']['maximumFoldWithinLibraryR2Loss'])))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'metadata', 'design', 'protocol', 'scores', 'program-reference', 'out']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--library-mean-only', action='store_true')
    args = parser.parse_args(); assert not args.out.exists()
    began = time.monotonic(); protocol = json.loads(args.protocol.read_text())
    inputs = {k: getattr(args, k) for k in ['source', 'metadata', 'design', 'protocol', 'scores']}
    inputs['programReference'] = args.program_reference/'reference.json'
    reference = json.loads(inputs['programReference'].read_text())
    assert reference['status'] == 'passed' and reference['cells'] == protocol['cells']
    inputs['programScores'] = args.program_reference/'scores.npy'
    if args.baseline:
        inputs['baseline'] = args.baseline
    bindings = {k: dict(bytes=p.stat().st_size, SHA256=sha(p)) for k, p in inputs.items()}
    for name in ['source', 'metadata', 'design']:
        assert bindings[name]['SHA256'] == protocol[name+'SHA256']
    assert reference['bindings']['source'] == bindings['source']['SHA256']
    assert reference['bindings']['protocol'] == bindings['protocol']['SHA256']
    assert reference['bindings']['definitions'] == protocol['definitionsSHA256']
    assert [v['definition']['id'] for v in reference['programs']] == protocol['programs']
    assert reference['payloads']['scores.npy'] == bindings['programScores']
    targets = np.load(inputs['programScores'], mmap_mode='r', allow_pickle=False)
    assert targets.shape == (protocol['cells'], len(protocol['programs'])) and targets.dtype == np.dtype('<f8')
    design = json.loads(args.design.read_text()); samples, ids = validate_design(design)
    assert len(samples) == protocol['libraries']
    codes = source_codes(args.source, args.metadata, samples, protocol['cells'], protocol['rowBatch'])
    stats = moments(args.scores, targets, codes, len(samples), protocol['components'], protocol['rowBatch'])
    result = evaluate(stats, codes, ids, design, protocol, args.library_mean_only)
    result.update(status='measured', bindings=bindings, protocol=protocol, libraryMeanOnly=args.library_mean_only,
                  evaluatorSHA256=sha(Path(__file__)), dependencySHA256=sha(Path(__file__).with_name('check_integration_response.py')),
                  scope=protocol['scope'])
    if args.baseline:
        baseline = json.loads(args.baseline.read_text())
        assert baseline['status'] == 'measured' and not baseline['libraryMeanOnly'] and 'comparisons' not in baseline
        assert baseline['evaluatorSHA256'] == result['evaluatorSHA256'] and baseline['dependencySHA256'] == result['dependencySHA256']
        for k in ['source', 'metadata', 'design', 'protocol', 'programReference', 'programScores']:
            assert baseline['bindings'][k] == bindings[k], ('Baseline identity', k)
        assert baseline['bindings']['scores']['SHA256'] == protocol['baselineScoresSHA256']
        result['comparisons'] = compare(baseline, result, protocol)
        result['allProgramGradientGatesPassed'] = all(all(c['gates'].values()) for c in result['comparisons'])
    else:
        assert not args.library_mean_only and bindings['scores']['SHA256'] == protocol['baselineScoresSHA256']
    assert all(p.stat().st_size == bindings[k]['bytes'] and sha(p) == bindings[k]['SHA256'] for k, p in inputs.items()), 'Input changed'
    result['seconds'] = time.monotonic()-began
    result['maximumResidentBytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if platform.system() == 'Darwin' else 1024)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as f:
        json.dump(result, f, sort_keys=True, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'cells', 'unavailableEmptyCellTargets', 'seconds', 'maximumResidentBytes']}))


if __name__ == '__main__':
    main()
