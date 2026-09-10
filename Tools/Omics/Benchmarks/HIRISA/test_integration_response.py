"""Independent dense/sklearn checks for the streamed response evaluator."""
import json
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

import check_integration_response as check


def fixture():
    rng = np.random.default_rng(91)
    samples, pairs, values, codes = [], [], [], []
    for donor in range(4):
        pair = dict(donor=str(donor), batch='batch'+str(donor), pool='p')
        for condition in ['none', 'drug']:
            k = len(samples); name = 'library'+str(k)
            samples.append({'accession': name, 'subject id': [str(donor)],
                            'cell type': ['enriched'], 'treatment': [condition],
                            'batch id': [pair['batch']], 'pool id': ['p'],
                            'batch pool': [pair['batch']+'-p']})
            pair['control' if condition == 'none' else 'treated'] = name
            n = 31+k
            x = rng.normal(size=(n, 3))
            x[:, 0] += 4*(condition == 'drug')
            x[:, 1] += donor/10
            values.append(x); codes.extend([k]*n)
        pairs.append(pair)
    # Shuffle all rows to ensure library boundaries and tile boundaries differ.
    x = np.concatenate(values); codes = np.asarray(codes, dtype=np.uint16)
    permutation = rng.permutation(len(x))
    return x[permutation], codes[permutation], dict(assignmentUsesExpression=False, samples=samples,
        comparisons=[dict(population='enriched', treatment='drug', pairs=pairs)])


def write_matrix(path, x):
    a = np.empty(x.shape, dtype=check.RECORD)
    a['row'] = np.arange(len(x))[:, None]; a['column'] = np.arange(x.shape[1]); a['value'] = x
    path.write_bytes(a.tobytes())


class ResponseChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
        self.x, self.codes, self.design = fixture()
        self.samples, self.ids = check.validate_design(self.design)
        self.path = self.root/'scores.bin'; write_matrix(self.path, self.x)

    def tearDown(self):
        self.tmp.cleanup()

    def evaluate(self, erasure=False):
        return check.evaluate(self.path, self.codes, self.samples, self.ids, self.design, 3, 17, 1., erasure)

    def test_dense_weighted_sklearn_all_folds_and_confusions(self):
        actual = self.evaluate()
        for f in actual['folds']:
            train_ids = [self.ids.index(v) for v in f['trainingLibraries']]
            mask = np.isin(self.codes, train_ids); x = self.x[mask]; cc = self.codes[mask]
            weights = np.array([1/(len(train_ids)*(cc == c).sum()) for c in cc])
            target = np.where(cc % 2 == 0, -1., 1.)
            scale = StandardScaler().fit(x, sample_weight=weights)
            model = Ridge(alpha=1., fit_intercept=True, solver='cholesky').fit(scale.transform(x), target, sample_weight=weights)
            np.testing.assert_allclose(f['weight'], model.coef_/scale.scale_, atol=2e-14, rtol=2e-13)
            np.testing.assert_allclose(f['intercept'], model.intercept_-scale.mean_@(model.coef_/scale.scale_), atol=2e-14, rtol=2e-13)
            for slot, role in enumerate(['control', 'treated']):
                held = self.x[self.codes == self.ids.index(f[role])]
                prediction = model.predict(scale.transform(held)) >= 0
                self.assertEqual(f['correct'][slot], int((prediction == slot).sum()))
        self.assertEqual(actual['classificationCells'], len(self.x))
        self.assertEqual(actual['cellsOutsideMatchedClassifiers'], 0)

    def test_dense_moments_all_rows_and_scatter(self):
        size, means, second = check.moments(self.path, self.codes, 8, 3, 17)
        for k in range(8):
            x = self.x[self.codes == k]
            self.assertEqual(size[k], len(x))
            np.testing.assert_allclose(means[k], x.mean(axis=0), atol=1e-14)
            np.testing.assert_allclose(second[k], x.T@x/len(x), atol=1e-14)
        result = self.evaluate()
        for s in result['conditionalDonorScatter']:
            ks = [k for k, sample in enumerate(self.samples) if sample['treatment'] == [s['condition']]]
            center = means[ks].mean(axis=0)
            expected = np.mean(np.sum((means[ks]-center)**2, axis=1))
            np.testing.assert_allclose(s['betweenDonorSquaredDistance'], expected, atol=1e-14)

    def test_identity_and_mean_erasure_controls(self):
        baseline = self.evaluate(); erased = self.evaluate(True)
        margins = dict(minimumResponseNorm=1e-10, maximumMeanAccuracyLoss=.05,
                       maximumFoldAccuracyLoss=.1, maximumMeanResponseDrift=.25)
        identity = check.compare(baseline, baseline, self.design, self.ids, margins)[0]
        self.assertTrue(all(identity['gates'].values()))
        control = check.compare(baseline, erased, self.design, self.ids, margins)[0]
        self.assertFalse(control['gates']['meanResponsePreserved'])
        self.assertFalse(control['gates']['meanClassificationPreserved'])
        self.assertEqual(control['meanRelativeResponseDrift'], 1.)
        self.assertTrue(all(f['balancedAccuracy'] == .5 for f in erased['folds']))

    def test_negligible_response_is_unavailable(self):
        baseline = self.evaluate(); baseline['libraryMeans'] = np.zeros((8, 3)).tolist()
        margins = dict(minimumResponseNorm=1e-10, maximumMeanAccuracyLoss=.05,
                       maximumFoldAccuracyLoss=.1, maximumMeanResponseDrift=.25)
        result = check.compare(baseline, baseline, self.design, self.ids, margins)[0]
        self.assertFalse(result['gates']['completeResponsePairs'])
        self.assertFalse(result['gates']['meanResponsePreserved'])

    def test_matrix_corruptions_rejected(self):
        original = self.path.read_bytes()
        for field, value in [('row', 9), ('column', 9), ('value', np.nan), ('value', np.inf)]:
            with self.subTest(field=field, value=value):
                a = np.frombuffer(original, dtype=check.RECORD).copy(); a[field][19] = value
                self.path.write_bytes(a.tobytes())
                with self.assertRaises(AssertionError): list(check.blocks(self.path, len(self.x), 3, 17))
        self.path.write_bytes(original[:-1])
        with self.assertRaises(AssertionError): list(check.blocks(self.path, len(self.x), 3, 17))
        self.path.write_bytes(original+b'0')
        with self.assertRaises(AssertionError): list(check.blocks(self.path, len(self.x), 3, 17))

    def test_metadata_source_identity_and_rejections(self):
        hpath = self.root/'source.h5ad'; meta = self.root/'metadata.json'
        cells = [dict(sampleID=self.ids[k], barcode='barcode'+str(i)) for i, k in enumerate(self.codes)]
        native = [dict(id=s['accession'], donorID=s['subject id'][0], condition=s['treatment'][0], batchID=s['batch pool'][0]) for s in self.samples]
        meta.write_text(json.dumps(dict(samples=native, cells=cells)))
        with h5py.File(hpath, 'w') as h:
            h.create_group('X').attrs['shape'] = [len(self.x), 3]
            obs = h.create_group('obs'); dtype = h5py.string_dtype()
            obs.create_dataset('_index', data=[c['barcode'] for c in cells], dtype=dtype)
            obs.create_dataset('geo_accession', data=[c['sampleID'] for c in cells], dtype=dtype)
            for field, key in [('geo_donor', 'subject id'), ('geo_treatment', 'treatment'), ('geo_enrichment', 'cell type'), ('geo_batch_pool', 'batch pool')]:
                obs.create_dataset(field, data=[self.samples[k][key][0] for k in self.codes], dtype=dtype)
        np.testing.assert_array_equal(check.source_codes(hpath, meta, self.samples, len(self.x), 17), self.codes)
        for field, value in [('barcode', 'wrong'), ('sampleID', 'wrong')]:
            changed = json.loads(meta.read_text()); changed['cells'][18][field] = value
            bad = self.root/'bad.json'; bad.write_text(json.dumps(changed))
            with self.assertRaises(AssertionError): check.source_codes(hpath, bad, self.samples, len(self.x), 17)
        with h5py.File(hpath, 'r+') as h: h['obs/geo_treatment'][18] = 'wrong'
        with self.assertRaises(AssertionError): check.source_codes(hpath, meta, self.samples, len(self.x), 17)

    def test_mismatched_pair_and_duplicate_donor_rejected(self):
        self.design['comparisons'][0]['pairs'][0]['batch'] = 'wrong'
        with self.assertRaises(AssertionError): check.validate_design(self.design)
        _, _, design = fixture(); design['comparisons'][0]['pairs'][1]['donor'] = '0'
        with self.assertRaises(AssertionError): check.validate_design(design)

    def test_baseline_binding_and_chained_result_rejected(self):
        bindings = {key: dict(bytes=1, SHA256=key) for key in ['source', 'metadata', 'design', 'protocol', 'scores']}
        baseline = dict(status='measured', libraryMeanErasure=False, bindings=bindings, evaluatorSHA256='code')
        protocol = dict(baselineScoresSHA256='scores')
        check.validate_baseline(baseline, bindings, 'code', protocol)
        for key in bindings:
            bad = json.loads(json.dumps(baseline)); bad['bindings'][key]['SHA256'] = 'wrong'
            with self.assertRaises(AssertionError): check.validate_baseline(bad, bindings, 'code', protocol)
        baseline['comparisons'] = []
        with self.assertRaises(AssertionError): check.validate_baseline(baseline, bindings, 'code', protocol)


if __name__ == '__main__':
    unittest.main(verbosity=2)
