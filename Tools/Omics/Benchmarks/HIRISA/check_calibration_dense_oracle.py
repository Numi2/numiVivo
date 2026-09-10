#!/usr/bin/env python3
"""Check post-result calibrated fits using centered per-cell SVD and residuals.

This independent oracle holds only the N-by-PC latent matrix densely. It does
not allocate cells by genes, and is not the streamed production evaluator.
"""
import argparse
import json
from pathlib import Path
import platform
import resource
import time
import numpy as np
from check_integration_response import blocks, sha, source_codes, validate_design


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'metadata', 'design', 'protocol', 'scores', 'program-reference', 'result', 'out', 'calibration-protocol']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args(); assert not args.out.exists()
    began = time.monotonic(); result = json.loads(args.result.read_text())
    assert sha(args.calibration_protocol)==result['calibrationProtocolSHA256']
    calibration=json.loads(args.calibration_protocol.read_text())
    protocol = json.loads(args.protocol.read_text()); design = json.loads(args.design.read_text())
    assert calibration['ridge']==protocol['ridge'] and calibration['margins']==protocol['margins']
    paths = dict(source=args.source, metadata=args.metadata, design=args.design, protocol=args.protocol,
                 scores=args.scores, programReference=args.program_reference/'reference.json',
                 programScores=args.program_reference/'scores.npy')
    for name, path in paths.items():
        assert result['bindings'][name] == dict(bytes=path.stat().st_size, SHA256=sha(path))
    result_sha = sha(args.result)
    assert not result['libraryMeanOnly'] and result['status'] == 'measured'
    samples, ids = validate_design(design)
    codes = source_codes(args.source, args.metadata, samples, protocol['cells'], protocol['rowBatch'])
    x = np.empty((protocol['cells'], protocol['components']))
    for first, last, values in blocks(args.scores, len(codes), x.shape[1], protocol['rowBatch']):
        x[first:last] = values
    y = np.load(args.program_reference/'scores.npy', mmap_mode='r', allow_pickle=False)
    assert np.isfinite(y).all(), 'This full-real oracle requires the recorded no-empty-cell cohort'
    index = {name: i for i, name in enumerate(ids)}
    maxima = dict(weight=0., intercept=0., withinLibraryR2=0., mse=0.)
    for fold in result['folds']:
        assert fold['status'] == 'measured'
        selected = [index[name] for name in fold['trainingLibraries']]
        train = np.isin(codes, selected)
        rows = np.flatnonzero(train)
        sizes = np.bincount(codes, minlength=len(ids))
        weights = 1/(len(selected)*sizes[codes[rows]])
        a = x[rows]; b = np.asarray(y[rows])
        center = np.average(a, axis=0, weights=weights)
        target_mean = np.average(b, axis=0, weights=weights)
        centered_x=a.copy();centered_y=b.copy()
        for k in selected:
            mask=codes[rows]==k
            centered_x[mask]-=a[mask].mean(axis=0)
            centered_y[mask]-=b[mask].mean(axis=0)
        scale = np.sqrt(np.average(centered_x**2, axis=0, weights=weights))
        scale[scale < 1e-12] = 1
        # Solve the augmented row system directly with SVD, without sufficient
        # statistics or the evaluator's normal equations.
        augmented = np.vstack([centered_x/scale*np.sqrt(weights[:, None]),
                               np.sqrt(protocol['ridge'])*np.eye(a.shape[1])])
        target = np.vstack([centered_y*np.sqrt(weights[:, None]), np.zeros((a.shape[1], b.shape[1]))])
        coefficients, _, rank, _ = np.linalg.lstsq(augmented, target, rcond=None)
        assert rank == a.shape[1]
        w = coefficients/scale[:, None]; intercept = target_mean-center@w
        for key, actual in [('weight', w), ('intercept', intercept)]:
            expected = np.asarray(fold[key]); error = float(np.max(np.abs(actual-expected)))
            maxima[key] = max(maxima[key], error)
            np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-10)
        for role, library in enumerate(fold['evaluationLibraries']):
            test = codes == index[library]; actual = np.asarray(y[test]); prediction = x[test]@w+intercept
            variance = actual.var(axis=0)
            centered_error = prediction-prediction.mean(axis=0)-(actual-actual.mean(axis=0))
            mse = ((prediction-actual)**2).mean(axis=0)
            for j, measured in enumerate(fold['measurements']):
                expected = measured['libraries'][role]
                if variance[j] <= protocol['margins']['minimumProgramVariance']:
                    assert expected['withinLibraryR2'] is None
                else:
                    skill = 1-float((centered_error[:, j]**2).mean()/variance[j])
                    error = abs(skill-expected['withinLibraryR2']); maxima['withinLibraryR2'] = max(maxima['withinLibraryR2'], error)
                    np.testing.assert_allclose(skill, expected['withinLibraryR2'], rtol=1e-9, atol=1e-9)
                maxima['mse'] = max(maxima['mse'], abs(float(mse[j])-expected['mse']))
                np.testing.assert_allclose(mse[j], expected['mse'], rtol=1e-9, atol=1e-10)
        del a, b, augmented, target
    assert sha(args.calibration_protocol)==result['calibrationProtocolSHA256']
    assert sha(args.result) == result_sha
    for name, path in paths.items():
        assert result['bindings'][name] == dict(bytes=path.stat().st_size, SHA256=sha(path))
    report = dict(status='passed', cells=len(codes), folds=len(result['folds']), programs=len(protocol['programs']),
                  maximumErrors=maxima, resultSHA256=result_sha, checkerSHA256=sha(Path(__file__)),
                  seconds=time.monotonic()-began,
                  maximumResidentBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if platform.system() == 'Darwin' else 1024),
                  scope='Every original held-out fold and program checked with full-cell within-library-centered weighted SVD and direct held-out residuals. Verifies the reported arithmetic, including negative skill and failed preservation margins; not biological calibration.')
    with args.out.open('x') as f:
        json.dump(report, f, indent=2, sort_keys=True); f.write('\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
