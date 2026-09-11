"""Independent numerical reconstruction and frozen external outcome comparison."""
from pathlib import Path
import argparse, gzip, hashlib, importlib.metadata, json
import anndata as ad
import numpy as np
from scipy import sparse
from scipy.stats import t
from sklearn.kernel_ridge import KernelRidge

root = Path(__file__).resolve().parent
src = root / 'prediction-inputs'
native = root / 'prediction-execution'
parser = argparse.ArgumentParser()
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(exist_ok=False)

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        while b := f.read(1048576): h.update(b)
    return h.hexdigest()

def write(name, d):
    (args.out / name).write_text(json.dumps(d, indent=2, sort_keys=True, allow_nan=False) + '\n')

freeze = json.loads((native / 'prediction-freeze.json').read_bytes())
assert not freeze['scoringStarted'] and freeze['donorPredictions'] == 124
assert freeze['inputFreezeSHA256'] == sha(src / 'input-freeze.json')
input_freeze = json.loads((src / 'input-freeze.json').read_bytes())
for base, records in [(src, input_freeze['files']), (native, freeze['files'])]:
    for p, digest in records.items(): assert sha(base / p) == digest, p
assert input_freeze['protocolSHA256'] == sha(root / 'PROTOCOL.md')
index = json.loads((native / 'retention.json').read_bytes())
for digest, info in index['objects'].items():
    p = native / 'objects' / (digest + '.gz')
    assert sha(p) == info['compressedSHA256'] and p.stat().st_size == info['compressedBytes']
    h = hashlib.sha256(); size = 0
    with gzip.open(p, 'rb') as f:
        while b := f.read(1048576): h.update(b); size += len(b)
    assert h.hexdigest() == digest and size == info['bytes']
assert {m['SHA256'] for b in index['bundles'].values() for m in b.values()} == set(index['objects'])

def get(bundle, name):
    item = index['bundles'][bundle][name]
    return json.loads(gzip.decompress((native / 'objects' / (item['SHA256'] + '.gz')).read_bytes()))

panel = json.loads((src / 'panel.json').read_bytes())
cohort = json.loads((src / 'cohort.json').read_bytes())
chunks = json.loads((src / 'chunks.json').read_bytes())
handoff = root / 'handoff'
assert sha(handoff / 'input-freeze.json') == input_freeze['originalHandoffFreezeSHA256']
for p, h in json.loads((handoff / 'input-freeze.json').read_bytes())['files'].items(): assert sha(handoff / p) == h
groups = json.loads((handoff / 'reference-groups.json').read_bytes())
counts = sparse.load_npz(handoff / 'reference-pseudobulk.npz').toarray().astype(np.uint64)
source_native = root / 'native-handoff-qualified'
assert sha(source_native / 'execution.json') == input_freeze['nativeAggregationExecutionSHA256']
assert sha(source_native / 'aggregate/report.json') == json.loads((source_native / 'execution.json').read_bytes())['reportSHA256']
bulk = json.loads((source_native / 'aggregate/report.json').read_bytes())['pseudobulk']
lookup = {(g['expID'], g['conditionCode']): i for i, g in enumerate(groups)}
bm = bulk['matrix']
native_counts = sparse.csr_matrix((bm['counts'], bm['featureIndices'], bm['rowOffsets']),
    shape=(bm['cellCount'], bm['featureCount'])).toarray().astype(np.uint64)
for row, g in enumerate(bulk['groups']):
    key = (g['donorID'].removeprefix('GSE181897:exp_id:'), g['condition'].removeprefix('source-code:'))
    np.testing.assert_array_equal(native_counts[row], counts[lookup[key]])
source_order = [bulk['featureIDs'].index(x.removeprefix('symbol|')) for x in panel]
source_logs = np.log1p(counts.astype(np.float64) / counts.sum(axis=1, keepdims=True, dtype=np.uint64) * 1e6)
scores, intervals, comparisons, model_checks = [], [], [], []
baselines = ['noChange', 'meanResponse', 'medianResponse', 'contextRidge']

for origin in ['Kang', 'HIRISA']:
    model = get('model-' + origin, 'model.json')
    assert model['featureIDs'] == panel
    training = ad.read_h5ad(src / origin / 'training.h5ad')
    fit = json.loads((src / origin / 'training.json').read_bytes())
    samples = {s['id']:s for s in fit['mapping']['samples']}
    rows = [samples[s] for s in training.obs['sample']]
    all_counts = training.X.toarray().astype(np.uint64)
    all_logs = np.log1p(all_counts.astype(float) / all_counts.sum(axis=1, keepdims=True, dtype=np.uint64) * 1e6)
    order = [training.var_names.get_loc(x) for x in panel]
    donors = sorted({s['donorID'] for s in rows})
    assert donors == model['trainingDonors'] and len(donors) == (8 if origin == 'Kang' else 5)
    pairs = []
    for d in donors:
        c = [i for i, s in enumerate(rows) if s['donorID'] == d and s['condition'] == 'control']
        y = [i for i, s in enumerate(rows) if s['donorID'] == d and s['condition'] == 'IFNB']
        assert len(c) == len(y) == 1
        pairs.append((c[0], y[0]))
    c, y = np.array(pairs).T
    controls = all_logs[c][:, order]
    response = all_logs[y][:, order] - controls
    pc = all_counts[c][:, order] + all_counts[y][:, order]
    selected = np.flatnonzero((pc.sum(axis=0, dtype=np.uint64) >= 10) & ((pc > 0).sum(axis=0) >= 2))
    assert selected.tolist() == model['selectedFeatureIndices'] and len(selected) > 0
    center = controls[:, selected].mean(axis=0)
    scale = controls[:, selected].std(axis=0)
    constant = np.all(controls[:, selected] == controls[0, selected], axis=0)
    center[constant] = controls[0, selected][constant]
    scale[constant | (scale == 0)] = 1
    contexts = (controls[:, selected] - center) / scale / np.sqrt(len(selected))
    mean, median = response.mean(axis=0), np.median(response, axis=0)
    kernel = contexts @ contexts.T
    dual = np.linalg.solve(kernel + np.eye(len(donors)), response - mean)
    estimator = KernelRidge(alpha=1, kernel='precomputed').fit(kernel, response - mean)
    np.testing.assert_allclose(estimator.dual_coef_, dual, rtol=1e-8, atol=1e-9)
    max_model = 0.
    for name, expected in [('contextCenters',center), ('contextScales',scale), ('contexts',contexts),
        ('meanResponse',mean), ('medianResponse',median), ('dualCoefficients',dual)]:
        actual = np.asarray(model[name])
        np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-9)
        max_model = max(max_model, float(np.max(np.abs(actual-expected))))
    variance = response.var(axis=0, ddof=1)
    available = (~np.all(response == response[0], axis=0)) & (variance > 0)
    native_variance = np.asarray([np.nan if x is None else x for x in model['donorResponseVariances']])
    assert np.array_equal(np.isfinite(native_variance), available)
    np.testing.assert_allclose(native_variance[available], variance[available], atol=1e-11, rtol=1e-9)
    critical = float(t.ppf(.975, len(donors)-1))
    half = critical * np.sqrt(variance * (1 + 1/len(donors)))
    model_checks.append(dict(origin=origin, donors=len(donors), selectedContextFeatures=len(selected), maximumDifference=max_model))
    seen = []
    for chunk in chunks:
        report = get('prediction-' + origin + '-' + chunk['id'], 'report.json')
        assert report['featureIDs'] == panel and len(report['predictions']) == len(chunk['donors'])
        expected_donors = {'GSE181897:exp_id:' + d for d in chunk['donors']}
        assert {p['group']['donorID'] for p in report['predictions']} == expected_donors
        for pred in report['predictions']:
            donor = pred['group']['donorID'].removeprefix('GSE181897:exp_id:')
            seen.append(donor)
            ctrl_row, truth_row = lookup[donor, 'C'], lookup[donor, 'B']
            control = source_logs[ctrl_row, source_order]
            truth = source_logs[truth_row, source_order]
            np.testing.assert_allclose(pred['control'], control, rtol=1e-10, atol=1e-10)
            assert pred['libraryCounts'] == int(counts[ctrl_row].sum(dtype=np.uint64))
            qc = (control[selected] - center) / scale / np.sqrt(len(selected))
            ridge = estimator.predict((qc @ contexts.T)[None,:])[0] + mean
            refs = dict(noChange=np.zeros(len(panel)), meanResponse=mean, medianResponse=median, contextRidge=ridge)
            assert [e['baseline'] for e in pred['estimates']] == baselines
            for e in pred['estimates']:
                baseline = e['baseline']; raw = control + refs[baseline]; expected = np.maximum(0, raw)
                for name, value in [('unclippedResponse',refs[baseline]),('predictedTreated',expected),('predictedResponse',expected-control)]:
                    np.testing.assert_allclose(e[name], value, rtol=1e-8, atol=1e-8)
                np.testing.assert_allclose(e['impliedCPMSum'], np.expm1(expected).sum(), rtol=1e-9, atol=1e-5)
                actual = np.asarray(e['predictedTreated']); error = actual - truth
                comparisons.append(dict(origin=origin, donor=donor, baseline=baseline, maximumDifference=float(np.max(np.abs(actual-expected)))))
                scores.append(dict(origin=origin, donor=donor, baseline=baseline, responseRMSE=float(np.sqrt(np.mean(error**2))),
                    responseMAE=float(np.mean(np.abs(error))), clippedFeatures=int(np.count_nonzero(raw < 0)),
                    controlCells=groups[ctrl_row]['cells'], treatedCells=groups[truth_row]['cells']))
            interval = pred['meanResponsePredictiveInterval']
            assert interval['nominalCoverage'] == .95 and interval['trainingDonors'] == len(donors) and interval['degreesOfFreedom'] == len(donors)-1
            assert interval['unavailableFeatureIndices'] == np.flatnonzero(~available).tolist()
            np.testing.assert_allclose(interval['studentCriticalValue'], critical, atol=1e-11, rtol=1e-11)
            refs_bounds = dict(unclippedResponseLower=mean-half, unclippedResponseUpper=mean+half,
                predictedTreatedLower=np.maximum(0, control+mean-half), predictedTreatedUpper=np.maximum(0, control+mean+half))
            bounds = {}; max_bound = 0.
            for name, expected in refs_bounds.items():
                actual = np.asarray([np.nan if x is None else x for x in interval[name]])
                assert np.array_equal(np.isfinite(actual), available)
                np.testing.assert_allclose(actual[available], expected[available], rtol=1e-9, atol=1e-9)
                bounds[name] = actual
                max_bound = max(max_bound, float(np.max(np.abs(actual[available]-expected[available]))))
            for space, observed in [('unclippedResponse', truth-control), ('predictedTreated',truth)]:
                lo, hi, obs = bounds[space+'Lower'][available], bounds[space+'Upper'][available], observed[available]
                intervals.append(dict(origin=origin, donor=donor, space=space, totalFeatures=len(panel),
                    availableFeatures=int(available.sum()), coveredFeatures=int(np.count_nonzero((obs>=lo)&(obs<=hi))),
                    belowFeatures=int(np.count_nonzero(obs<lo)), aboveFeatures=int(np.count_nonzero(obs>hi)),
                    coverage=float(np.mean((obs>=lo)&(obs<=hi))), meanWidth=float(np.mean(hi-lo)), maximumBoundDifference=max_bound))
    assert sorted(seen, key=int) == cohort['eligibleDonors'] and len(seen) == len(set(seen)) == 62

summaries = []
for origin in ['Kang','HIRISA']:
    means = {b:float(np.mean([s['responseRMSE'] for s in scores if s['origin']==origin and s['baseline']==b])) for b in baselines}
    by_donor = {d:{s['baseline']:s['responseRMSE'] for s in scores if s['origin']==origin and s['donor']==d} for d in cohort['eligibleDonors']}
    gain = 100 * (1-means['meanResponse']/means['noChange'])
    summaries.append(dict(origin=origin, donors=62, meanRMSE=means, meanGainPercent=gain, requiredMeanGainPercent=5,
        primaryPass=gain>=5, ridgeSecondaryPass=means['contextRidge']<min(means['noChange'],means['meanResponse']),
        meanWorseNoChange=sum(x['meanResponse']>x['noChange'] for x in by_donor.values()),
        ridgeWorseNoChange=sum(x['contextRidge']>x['noChange'] for x in by_donor.values()),
        ridgeWorseMean=sum(x['contextRidge']>x['meanResponse'] for x in by_donor.values()),
        intervalSummary={space:dict(coverageMean=float(np.mean([x['coverage'] for x in intervals if x['origin']==origin and x['space']==space])),
            coverageRange=[float(fn([x['coverage'] for x in intervals if x['origin']==origin and x['space']==space])) for fn in [min,max]],
            meanWidth=float(np.mean([x['meanWidth'] for x in intervals if x['origin']==origin and x['space']==space])))
            for space in ['unclippedResponse','predictedTreated']}))
assert len(scores)==496 and len(intervals)==248 and len(comparisons)==496
write('scores.json', scores); write('intervals.json', intervals); write('comparisons.json', comparisons)
write('model-checks.json', model_checks); write('summary.json', summaries)
write('checks.json', dict(status='passed-independent-reconstruction-and-complete-scoring', donorPredictions=124,
    vectors=496, intervalAssessments=248, features=11800, allSourceNativeAggregateRowsCompared=379,
    inputFreezeSHA256=sha(src/'input-freeze.json'), predictionFreezeSHA256=sha(native/'prediction-freeze.json'),
    scorerSHA256=sha(__file__), packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','scikit-learn','anndata']},
    maximumPredictionDifference=max(x['maximumDifference'] for x in comparisons),
    maximumIntervalDifference=max(x['maximumBoundDifference'] for x in intervals)))
print(json.dumps(summaries, indent=2))
