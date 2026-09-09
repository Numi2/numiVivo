#!/usr/bin/env python3
"""Check complete replay and fold isolation using the real Norman reference."""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from combinations import load, sha, write
from unseen_targets import METHODS, fit_predict, select

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference', type=Path, required=True)
p.add_argument('--first', type=Path, required=True)
p.add_argument('--repeat', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
assert not a.output.exists()
first, repeat = load(a.first / 'predictions.npz'), load(a.repeat / 'predictions.npz')
assert first.keys() == repeat.keys()
for k in first:
    np.testing.assert_array_equal(first[k], repeat[k])
for file in ('folds.json', 'receipt.json'):
    assert (a.first / file).read_bytes() == (a.repeat / file).read_bytes()
folds = json.loads((a.first / 'folds.json').read_text())
assert len(folds) == 105 and int(first['supported'].sum()) == 102
assert [r['target'] for r in folds if r['descriptorStatus'] != 'supported'] == ['C19orf26', 'C3orf72', 'KIAA1804']
for i, r in enumerate(folds):
    assert r['target'] not in r['trainingTargets'] and len(r['trainingTargets']) == 104
    assert r['excludedOutcomeMutationInputsExact']
    for method in METHODS:
        v = first[method][i]
        if method in METHODS[2:] and not first['supported'][i]:
            assert np.all(np.isnan(v))
        else:
            assert np.all(np.isfinite(v) & (v >= 0))
ref = load(a.reference)
data = select(ref, 'AHR')
with threadpool_limits(limits=1):
    expected = fit_predict(data)
    changed = dict(ref, counts=ref['counts'].copy())
    names = ref['conditions'].tolist()
    changed['counts'][[names.index(n) for n in data['excluded']]] = 1
    actual = fit_predict(select(changed, 'AHR'))
    for method in expected[0]:
        np.testing.assert_array_equal(expected[0][method], actual[0][method])
    np.testing.assert_array_equal(expected[2], actual[2])
    # Source condition ordering cannot change the canonical training selection.
    reordered = dict(ref, conditions=ref['conditions'][::-1], counts=ref['counts'][::-1])
    other = select(reordered, 'AHR')
    for key in ('control', 'trainingCounts'):
        np.testing.assert_array_equal(data[key], other[key])
    rejected = []
    for name, value in [('negative', -1), ('fractional', 0.5), ('nonfinite', float('nan')), ('inexact', float(2**54))]:
        bad = copy.deepcopy(data)
        bad['trainingCounts'] = bad['trainingCounts'].astype(float)
        bad['trainingCounts'][0, 0] = value
        try:
            fit_predict(bad)
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError('accepted ' + name)
    bad = copy.deepcopy(data)
    bad['control'][:] = 0
    try:
        fit_predict(bad)
    except ValueError:
        rejected.append('empty-control')
    else:
        raise AssertionError('accepted empty control')
write(a.output, dict(status='passed', arraysReplayedExactly=list(first), folds=105,
    supportedTargets=102, excludedMutationInputFolds=105, excludedMutationPredictionReplayTarget='AHR',
    sourceConditionOrderExact=True, rejected=rejected, referenceSHA256=sha(a.reference),
    predictionReceiptSHA256=sha(a.first / 'receipt.json'), implementationSHA256=sha(Path(__file__))))
