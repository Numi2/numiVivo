#!/usr/bin/env python3
"""Compare native sparse reduction with Scanpy on the exact imported count axes.

No cell-by-gene dense array is constructed. Dense arrays contain only PCA scores,
loadings, or the small cross-product used for sign/rotation invariant comparison.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import anndata as ad
import numpy as np
import scanpy as sc
import pandas as pd
from fast_array_utils import stats
from scipy import sparse

p = argparse.ArgumentParser()
source = p.add_mutually_exclusive_group(required=True)
source.add_argument('--dataset', type=Path)
source.add_argument('--h5ad', type=Path)
p.add_argument('--report', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
report = json.loads(a.report.read_text())
r = report['reduction']
if a.dataset:
    d = json.loads(a.dataset.read_text())
    m = d['matrix']
    x = sparse.csr_matrix((np.asarray(m['counts'], dtype=np.float64), m['featureIndices'], m['rowOffsets']),
                          shape=(m['cellCount'], m['featureCount']))
    identities = [(c['sampleID'], c['barcode']) for c in d['cells']]
    feature_ids = [f['id'] for f in d['features']]
    target = report['plan']['normalizationTarget']
    source_path = a.dataset
else:
    import h5py
    plan = json.loads((a.report.parent/'plan.json').read_text())
    mapping = plan['mapping']
    path = mapping['matrixPath']
    with h5py.File(a.h5ad, 'r') as handle:
        assert isinstance(handle[path], h5py.Group), 'Reference requires a sparse H5AD count array'
    obj = ad.read_h5ad(a.h5ad)
    if path == 'X':
        x, var = obj.X, obj.var
    elif path == 'raw/X':
        x, var = obj.raw.X, obj.raw.var
    else:
        assert path.startswith('layers/')
        x, var = obj.layers[path[len('layers/'):]], obj.var
    assert sparse.issparse(x)
    x = x.astype(np.float64).tocsr()
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    barcodes = obj.obs[mapping['barcodeColumn']].astype(str).tolist() if mapping.get('barcodeColumn') else obj.obs_names.tolist()
    identities = list(zip(obj.obs[mapping['sampleColumn']].astype(str), barcodes))
    feature_ids = var[mapping['featureIDColumn']].astype(str).tolist() if mapping.get('featureIDColumn') else var.index.tolist()
    target = report['reductionStorage']['options']['normalizationTarget']
    source_path = a.h5ad
assert len(set(identities)) == len(identities)
identity = {cell: i for i, cell in enumerate(identities)}
rows = [identity[(c['sampleID'], c['barcode'])] for c in r['cells']]
b = ad.AnnData(x[rows].copy())
b.var_names = feature_ids
assert list(b.var_names) == [f['featureID'] for f in r['features']]
sc.pp.normalize_total(b, target_sum=target)
sc.pp.log1p(b)
o = r['options']
sc.pp.highly_variable_genes(b, flavor='seurat', n_top_genes=o['highlyVariableFeatures'], n_bins=o['meanBins'])
selected = np.flatnonzero(b.var['highly_variable'].to_numpy())
assert selected.tolist() == r['selectedFeatureIndices'], 'HVG membership mismatch'
errors = {}
linear = b.X.copy()
linear.data = np.expm1(linear.data)
reference_mean, reference_variance = stats.mean_var(linear, axis=0, correction=1)
for native, expected in [('meanNormalized', reference_mean), ('varianceNormalized', reference_variance)]:
    observed = np.array([f[native] for f in r['features']])
    np.testing.assert_allclose(observed, expected, rtol=1e-9, atol=1e-9)
    errors[native] = float(np.max(np.abs(observed-expected)))
reference_bins = pd.cut(b.var['means'], bins=o['meanBins']).cat.codes.to_numpy()
assert reference_bins.tolist() == [f['meanBin'] for f in r['features']], 'mean bin mismatch'

for native, reference in [('logMeanForBinning', 'means'), ('logDispersion', 'dispersions'), ('normalizedDispersion', 'dispersions_norm')]:
    lhs = np.array([np.nan if f.get(native) is None else f[native] for f in r['features']])
    rhs = b.var[reference].to_numpy(dtype=np.float64)
    assert np.array_equal(np.isnan(lhs), np.isnan(rhs)), (native, 'undefined mask')
    valid = np.isfinite(lhs) & np.isfinite(rhs)
    # Scanpy stores normalized dispersion as float32; moments remain float64.
    tolerance = 2e-6 if native == 'normalizedDispersion' else 1e-9
    np.testing.assert_allclose(lhs, rhs, rtol=tolerance, atol=tolerance, equal_nan=True)
    errors[native] = float(np.max(np.abs(lhs[valid] - rhs[valid])))
c = b[:, selected].copy()
sc.pp.pca(c, n_comps=o['components'], zero_center=True, svd_solver='arpack', dtype='float64', random_state=7)
assert sparse.issparse(c.X)
loading = np.asarray(r['loadings'])
scores = np.asarray(r['scores'])
reference_loading = np.asarray(c.varm['PCs'])
variance = np.asarray(r['explainedVariance'])
np.testing.assert_allclose(variance, c.uns['pca']['variance'], rtol=1e-7, atol=1e-9)
np.testing.assert_allclose(r['explainedVarianceRatio'], c.uns['pca']['variance_ratio'], rtol=1e-7, atol=1e-10)
# Orthogonal Procrustes permits sign and within-subspace rotation ambiguity.
u, singular, vt = np.linalg.svd(loading.T @ reference_loading)
rotation = u @ vt
score_error = float(np.linalg.norm(scores @ rotation - c.obsm['X_pca']) / np.linalg.norm(c.obsm['X_pca']))
assert score_error < 1e-5, score_error
assert float(np.min(singular)) > 1 - 1e-8, singular
np.testing.assert_allclose(scores.mean(axis=0), 0, atol=1e-10)
np.testing.assert_allclose(scores.T @ scores / (len(rows)-1), np.diag(variance), rtol=1e-7, atol=1e-8)
assert max(r['relativeResiduals']) <= o['relativeResidualTolerance']
result = dict(status='passed', qualification='Numerical agreement on these count axes; no biological or integration qualification',
    datasetSHA256=hashlib.sha256(source_path.read_bytes()).hexdigest(), reportSHA256=hashlib.sha256(a.report.read_bytes()).hexdigest(),
    cells=len(rows), features=x.shape[1], sourceNonzeros=x.nnz, selectedFeatures=len(selected), components=o['components'],
    featureMaximumAbsoluteErrors=errors, minimumSubspaceCosine=float(min(singular)), alignedRelativeScoreError=score_error,
    maximumNativeResidual=max(r['relativeResiduals']), maximumVarianceRelativeError=float(np.max(np.abs(variance/c.uns['pca']['variance']-1))),
    versions={name: importlib.metadata.version(name) for name in ['scanpy','anndata','numpy','scipy','pandas']})
if 'reductionStorage' in report:
    storage = report['reductionStorage']
    entries = int(x[rows][:, selected].nnz)
    visits = entries * (2 * min(o['maximumBasis'], len(selected)) + 3 * o['components'])
    assert storage['sourcePasses'] == 3
    assert storage['selectedEntries'] == entries and storage['cacheBytes'] == 16 * entries
    assert storage['entryVisits'] == visits <= storage['options']['maximumEntryVisits']
    assert storage['cacheBytes'] <= storage['options']['maximumCacheBytes']
    result['streamedStorage'] = storage
a.out.parent.mkdir(parents=True, exist_ok=True)
a.out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
print(json.dumps(result, indent=2))
