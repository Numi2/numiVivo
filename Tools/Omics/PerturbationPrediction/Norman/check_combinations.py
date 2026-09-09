#!/usr/bin/env python3
"""Check leakage boundaries and frozen predictions using the full real input."""
import argparse
import json
from pathlib import Path

import numpy as np

import combinations as owner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reference', 'inputs', 'predictions', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    assert owner.sha(args.reference) == owner.REFERENCE_SHA256
    reference = owner.load(args.reference)
    original, queries, truth = owner.select_inputs(reference)
    altered = {key: value.copy() for key, value in reference.items()}
    rows = np.array(['_' in name for name in reference['conditions']])
    altered['counts'][rows] = altered['counts'][rows, ::-1] + 17
    training, mutated_queries, mutated_truth = owner.select_inputs(altered)
    assert queries == mutated_queries
    assert not np.array_equal(truth['counts'], mutated_truth['counts'])
    saved_training = owner.load(args.inputs / 'training.npz')
    for key in original:
        np.testing.assert_array_equal(original[key], training[key])
        np.testing.assert_array_equal(original[key], saved_training[key])
    np.testing.assert_array_equal(owner.load(args.inputs / 'sealed-outcomes.npz')['counts'], truth['counts'])
    frozen = owner.load(args.predictions / 'model.npz')
    model, _ = owner.fit(training)
    for key in model:
        np.testing.assert_array_equal(model[key], frozen[key])
    # Canonical training-target order must also preserve the mean/shuffle model.
    reordered = {key: value.copy() for key, value in training.items()}
    reordered['targetIDs'] = reordered['targetIDs'][::-1]
    reordered['singleCounts'] = reordered['singleCounts'][::-1]
    reordered_model, _ = owner.fit(reordered)
    for key in model:
        np.testing.assert_array_equal(model[key], reordered_model[key])
    predicted = {method: owner.load(args.predictions / (method + '.npz')) for method in owner.METHODS}
    comparisons = 0
    for i, query in enumerate(queries):
        values, _, _ = owner.predict_one(model, query['targets'])
        for method, value in values.items():
            np.testing.assert_array_equal(value, predicted[method]['expression'][predicted[method]['queryRows'][i]])
            comparisons += 1
    rejected = []
    for name, query in [('unknown-target', ['unobserved_target', 'AHR']), ('repeated-target', ['AHR', 'AHR']),
                        ('one-target', ['AHR']), ('three-targets', ['AHR', 'ARID1A', 'BCL2L11'])]:
        try:
            owner.predict_one(model, query)
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError(name)
    for name, value in [('negative-count', -1.0), ('fractional-count', 0.5), ('nan-count', np.nan),
                        ('infinite-count', np.inf), ('outside-exact-integer-range', np.uint64(2**53 + 1))]:
        array = np.array([value, 1], dtype=np.uint64 if isinstance(value, np.uint64) else float)
        try:
            owner.counts(array)
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError(name)
    try:
        owner.counts(np.zeros(10))
    except ValueError:
        rejected.append('empty-library')
    else:
        raise AssertionError('empty-library')
    args.out.mkdir(parents=True, exist_ok=False)
    result = dict(status='passed', sealedPairOutcomesMutated=131, trainingArraysUnchanged=True,
                  modelArraysUnchanged=True, frozenPredictionArraysExact=comparisons,
                  trainingTargetOrderExact=True, rejected=rejected,
                  referenceSHA256=owner.sha(args.reference),
                  predictionReceiptSHA256=owner.sha(args.predictions / 'receipt.json'))
    owner.write(args.out / 'checks.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
