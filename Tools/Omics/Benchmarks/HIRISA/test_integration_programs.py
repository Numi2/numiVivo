#!/usr/bin/env python3
"""Independent dense and weighted-regression controls for program diagnostics."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy import sparse
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from check_integration_response import RECORD
from check_integration_programs import moments, fit, library_metrics, evaluate, compare
from prepare_program_reference import resolve, score_block


class Programs(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        rng = np.random.default_rng(20260910)
        codes = np.repeat(np.arange(8), [7, 9, 13, 8, 12, 5, 11, 6])
        rng.shuffle(codes)
        self.codes = codes
        self.x = rng.normal(size=(len(codes), 3))+codes[:, None]*.2
        self.x[:, 2] = 4
        self.y = self.x @ np.array([[1.2, -.8], [.4, 1.1], [0, 0]])+rng.normal(scale=.1, size=(len(codes), 2))
        records = np.empty(self.x.shape, dtype=RECORD)
        records['row'] = np.arange(len(codes))[:, None]
        records['column'] = np.arange(3)
        records['value'] = self.x
        self.path = self.root/'scores.bin'; records.tofile(self.path)
        self.stats = moments(self.path, self.y, codes, 8, 3, 5)
        self.ids = ['library'+str(i) for i in range(8)]
        self.design = dict(comparisons=[dict(population='declared', treatment='known', pairs=[
            dict(donor=str(i), control=self.ids[2*i], treated=self.ids[2*i+1]) for i in range(4)])])
        self.protocol = dict(programs=['first', 'second'], ridge=1., margins={
            'minimumProgramVariance': 1e-12, 'minimumErasureDetection': .05,
            'maximumMeanWithinLibraryR2Loss': .05, 'maximumFoldWithinLibraryR2Loss': .1})

    def tearDown(self):
        self.temporary.cleanup()

    def test_sparse_scores_and_empty_library(self):
        counts = np.array([[3, 0, 2, 0], [0, 0, 0, 0], [0, 4, 0, 1], [8, 2, 1, 0]], dtype=np.uint16)
        weights = np.array([[.5, 0], [0, 1], [-.5, 0], [0, 0]])
        csr = sparse.csr_matrix(counts)
        total, observed, detected, _ = score_block(csr.data, csr.indices, csr.indptr, weights, 10000)
        expected = np.zeros_like(counts, dtype=float); valid = counts.sum(axis=1) > 0
        expected[valid] = np.log1p(counts[valid]/counts[valid].sum(axis=1)[:, None]*10000)
        expected = expected @ weights; expected[~valid] = np.nan
        np.testing.assert_allclose(observed, expected, rtol=1e-14, atol=1e-14, equal_nan=True)
        np.testing.assert_array_equal(total, counts.sum(axis=1))
        np.testing.assert_array_equal(detected, (counts > 0).astype(int) @ (weights != 0).astype(int))
        with self.assertRaisesRegex(AssertionError, 'canonical'):
            score_block(np.array([1, 2], dtype=np.uint16), np.array([0, 0]), np.array([0, 2]), weights, 10000)

    def test_exact_symbol_resolution_and_coverage(self):
        definition = dict(featureNamespace='HUMAN_GENE_SYMBOL', members=[dict(featureID='A', weight=1), dict(featureID='missing', weight=1)], minimumWeightCoverage=.5)
        weights, resolved = resolve(['idB', 'idA'], ['B', 'A'], [definition])
        np.testing.assert_array_equal(weights[:, 0], [0, 1])
        self.assertEqual(resolved[0]['missingSymbols'], ['missing'])
        self.assertEqual(resolved[0]['weightCoverage'], .5)
        with self.assertRaisesRegex(AssertionError, 'Ambiguous'):
            resolve(['id1', 'id2'], ['A', 'A'], [definition])
        with self.assertRaises(AssertionError):
            resolve(['idB'], ['B'], [definition])

    def test_all_moments_against_dense(self):
        for k in range(8):
            x = self.x[self.codes == k]; y = self.y[self.codes == k]
            expected = dict(x=x.mean(axis=0), xx=x.T@x/len(x), y=y.mean(axis=0),
                            yy=(y*y).mean(axis=0), xy=x.T@y/len(x))
            for name, value in expected.items():
                np.testing.assert_allclose(self.stats[name][k], value, rtol=1e-13, atol=1e-13)
        for rows in [1, 13, 8192]:
            alternate = moments(self.path, self.y, self.codes, 8, 3, rows)
            for name, value in self.stats.items():
                np.testing.assert_allclose(alternate[name], value, rtol=1e-13, atol=1e-13)

    def test_every_fold_against_weighted_sklearn(self):
        for donor in range(4):
            selected = [i for i in range(8) if i//2 != donor]
            train = np.isin(self.codes, selected)
            sample_weight = np.array([1/(len(selected)*np.sum(self.codes == c)) for c in self.codes[train]])
            scaler = StandardScaler().fit(self.x[train], sample_weight=sample_weight)
            oracle = Ridge(alpha=1.).fit(scaler.transform(self.x[train]), self.y[train], sample_weight=sample_weight)
            weight, intercept, mean = fit(self.stats, selected, 1.)
            np.testing.assert_allclose(self.x@weight+intercept, oracle.predict(scaler.transform(self.x)), rtol=1e-12, atol=1e-12)
            for k in [2*donor, 2*donor+1]:
                chosen = self.codes == k; truth = self.y[chosen]; prediction = oracle.predict(scaler.transform(self.x[chosen]))
                measured = library_metrics(self.stats, k, weight, intercept, mean, 1e-12)
                centered_error = (prediction-prediction.mean(axis=0))-(truth-truth.mean(axis=0))
                expected_skill = 1-(centered_error**2).mean(axis=0)/truth.var(axis=0)
                expected_mse = ((prediction-truth)**2).mean(axis=0)
                for j in range(2):
                    self.assertAlmostEqual(measured[j]['withinLibraryR2'], expected_skill[j], places=11)
                    self.assertAlmostEqual(measured[j]['mse'], expected_mse[j], places=11)

    def test_identity_and_gradient_erasure(self):
        base = evaluate(self.stats, self.codes, self.ids, self.design, self.protocol)
        identity = evaluate(self.stats, self.codes, self.ids, self.design, self.protocol)
        self.assertEqual(base, identity)
        for value in compare(base, identity, self.protocol):
            self.assertTrue(all(value['gates'].values()))
        erased = evaluate(self.stats, self.codes, self.ids, self.design, self.protocol, True)
        for fold in erased['folds']:
            for measurement in fold['measurements']:
                self.assertEqual(measurement['withinLibraryR2'], 0.)
        for value in compare(base, erased, self.protocol):
            self.assertTrue(value['controlSensitive'])
            self.assertFalse(value['gates']['meanPreserved'])

    def test_missing_and_constant_targets_are_unavailable(self):
        values = self.y.copy(); values[0] = np.nan
        stats = moments(self.path, values, self.codes, 8, 3, 5)
        result = evaluate(stats, self.codes, self.ids, self.design, self.protocol)
        self.assertEqual(result['unavailableEmptyCellTargets'], 1)
        values[1, 0] = np.nan
        with self.assertRaisesRegex(AssertionError, 'Partially missing'):
            moments(self.path, values, self.codes, 8, 3, 5)
        stats = moments(self.path, np.ones_like(self.y), self.codes, 8, 3, 5)
        result = evaluate(stats, self.codes, self.ids, self.design, self.protocol)
        self.assertTrue(all(v['withinLibraryR2'] is None for f in result['folds'] for v in f['measurements']))
        self.assertTrue(all(not c['gates']['completeFolds'] for c in compare(result, result, self.protocol)))


if __name__ == '__main__':
    unittest.main(verbosity=2)
