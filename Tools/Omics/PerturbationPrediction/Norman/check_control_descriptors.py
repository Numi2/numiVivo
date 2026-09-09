#!/usr/bin/env python3
"""Real-data replay, coverage and exclusion checks for control descriptors."""
import argparse
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from combinations import load, sha, write
from control_descriptors import METHODS, fit_predict, select, verify_reference

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference', type=Path, required=True)
p.add_argument('--descriptors', type=Path, required=True)
p.add_argument('--repeat-descriptors', type=Path, required=True)
p.add_argument('--first', type=Path, required=True)
p.add_argument('--repeat', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
assert not a.out.exists()
for first, repeat, files in [(a.descriptors, a.repeat_descriptors, ['descriptors.npz', 'coverage.json', 'receipt.json']),
                             (a.first, a.repeat, ['predictions.npz', 'folds.json', 'receipt.json'])]:
    for f in files:
        assert sha(first / f) == sha(repeat / f), f
ref = verify_reference(a.reference)
d = load(a.descriptors / 'descriptors.npz')
pred = load(a.first / 'predictions.npz')
np.testing.assert_array_equal(pred['featureIDs'], ref['feature_ids'])
np.testing.assert_array_equal(d['controlBarcodes'], ref['barcodes'][ref['cell_conditions'] == 'control'])
assert len(d['landmarkIDs']) == 2000 and len(set(d['landmarkIDs'])) == 2000
assert np.all(np.isfinite(d['expression'])) and np.max(abs(d['expression'])) <= 1 + 1e-10
coverage = json.loads((a.descriptors / 'coverage.json').read_text())
assert len(coverage) == 105
assert [r['target'] for r in coverage if r['status'] == 'supported'] == d['targetIDs'].tolist()
features, symbols = ref['feature_ids'].tolist(), ref['gene_symbols'].tolist()
alltargets = pred['targetIDs'].tolist()
assert all(symbols[features.index(f)] not in alltargets for f in d['landmarkIDs'])
for i, target in enumerate(alltargets):
    assert bool(pred['supported'][i]) == (target in d['targetIDs'])
    for method in METHODS:
        if method in METHODS[2:] and not pred['supported'][i]:
            assert np.all(np.isnan(pred[method][i]))
        else:
            assert np.all(np.isfinite(pred[method][i]) & (pred[method][i] >= 0))
    data = select(ref, target)
    changed = dict(ref, counts=ref['counts'].copy())
    names = ref['conditions'].tolist()
    changed['counts'][[names.index(n) for n in data['excluded']]] = 1
    mutated = select(changed, target)
    for key in data:
        if isinstance(data[key], np.ndarray):
            np.testing.assert_array_equal(data[key], mutated[key])
        else:
            assert data[key] == mutated[key]
# Complete refits after excluded-outcome mutations at fixed early/middle/late targets.
with threadpool_limits(limits=1):
    for target in [alltargets[0], alltargets[len(alltargets)//2], alltargets[-1]]:
        data = select(ref, target)
        altered = dict(ref, counts=ref['counts'].copy())
        altered['counts'][[names.index(n) for n in data['excluded']]] = 1
        actual, baseline, top, info = fit_predict(select(altered, target), d)
        i = alltargets.index(target)
        for method, v in actual.items():
            np.testing.assert_array_equal(v, pred[method][i])
        np.testing.assert_array_equal(baseline, pred['baseline'][i])
        np.testing.assert_array_equal(top, pred['trainingTop1000'][i])
write(a.out, dict(status='passed', descriptorReplayByteExact=True, predictionReplayByteExact=True,
    foldMutationInputsExact=105, mutationRefitTargets=[alltargets[0], alltargets[len(alltargets)//2], alltargets[-1]],
    controls=11855, landmarks=2000, supportedTargets=int(pred['supported'].sum()),
    landmarkTargetsExcluded=True, unsupportedTargets=[r for r in coverage if r['status'] != 'supported'],
    predictionReceiptSHA256=sha(a.first / 'receipt.json'), descriptorReceiptSHA256=sha(a.descriptors / 'receipt.json'),
    implementationSHA256=sha(Path(__file__))))
