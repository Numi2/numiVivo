#!/usr/bin/env python3
"""Actual AnnData -> native Swift -> AnnData tests; not a biology benchmark."""
import argparse
import hashlib
import json
from importlib.metadata import version
from pathlib import Path
import shutil
import subprocess
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
checks = []

def run(name, source, plan, success=True):
    planpath = a.out / (name + '-plan.json')
    planpath.write_text(json.dumps(plan))
    result = subprocess.run([str(a.binary.resolve()), str(source), str(planpath), str(a.out / name)], capture_output=True, text=True)
    (a.out / (name + '.log')).write_text(result.stdout + result.stderr)
    assert result.returncode == (0 if success else 65), (name, result.returncode, result.stderr)
    checks.append(name)
    return a.out / name

samples = [dict(id=s, biologicalReplicateID=s, donorID='donor-' + s, condition='control' if s == 's1' else 'treated', batchID='batch', organism='NCBITaxon:9606') for s in ['s1', 's2']]
plan = dict(schemaVersion=1, id='interop', evidence='synthetic', sourceDescription='AnnData interoperability fixture, not biological measurements', countUnit='umiCount', matrixPath='layers/counts', samples=samples, sampleColumn='sample', groupColumn='group', featureNameColumn='name', mitochondrialFeatureIDs=['g1'])
counts = sparse.csr_matrix(np.array([[9007199254740993, 0, 3], [0, 4, 0], [5, 0, 0], [0, 0, 0]], dtype=np.uint64))
obs = pd.DataFrame({'sample': pd.Categorical(['s1','s1','s2','s2'], categories=['s2','s1'], ordered=True), 'group': pd.Categorical(['T','B',None,'T']), 'nullable': pd.array([1,None,2,3], dtype='Int64')}, index=pd.Index(['c1','c2','c3','c4'], dtype=object))
var = pd.DataFrame({'name': np.array(['α','same','same'], dtype=object)}, index=pd.Index(['g1','g2','g3'], dtype=object))
obj = ad.AnnData(X=sparse.csr_matrix(np.eye(4, 3, dtype=np.float32) * 0.5), obs=obs, var=var)
obj.layers['counts'] = counts
obj.obsm['X_pca'] = np.arange(8, dtype=np.float32).reshape(4,2)
obj.obsp['connectivities'] = sparse.eye(4, format='csr')
obj.uns['nested'] = {'labels': np.array(['kept','λ'], dtype=object), 'threshold': 0.25}
obj.raw = obj.copy()
for layout in ['csr','csc','dense']:
    obj.layers['counts'] = counts.toarray() if layout == 'dense' else counts.asformat(layout)
    source = a.out / (layout + '.h5ad')
    obj.write_h5ad(source)
    out = run(layout, source, plan)
    data = json.loads((out / 'dataset.json').read_text())
    assert data['matrix']['counts'] == [9007199254740993,3,4,5]
    assert data['matrix']['rowOffsets'] == [0,2,3,4,4]
    assert data['cells'][2].get('group') is None
    assert (out / 'original.h5ad').read_bytes() == source.read_bytes()
    restored = ad.read_h5ad(out / 'original.h5ad')
    pd.testing.assert_frame_equal(restored.obs, ad.read_h5ad(source).obs)
    assert restored.raw is not None and np.array_equal(restored.obsm['X_pca'], obj.obsm['X_pca'])
    native = ad.read_h5ad(out / 'native.h5ad')
    assert np.array_equal(native.X.toarray(), counts.toarray())
    assert native.obs['barcode'].tolist() == obj.obs_names.tolist()
    assert pd.isna(native.obs['group'].iloc[2])
    assert native.var['name'].tolist() == ['α','same','same']
    assert json.loads(native.uns['numivivo'])['samples'] == samples
    second = run(layout + '-native-reimport', out / 'native.h5ad', dict(plan, matrixPath='X', barcodeColumn='barcode'))
    assert json.loads((second / 'dataset.json').read_text()) == data

source = a.out / 'csr.h5ad'
# Noncanonical but legal SciPy CSR: duplicate entries, explicit zero and unsorted indices.
noncanonical = a.out / 'noncanonical.h5ad'
obj.layers['counts'] = sparse.csr_matrix((np.array([3, 2, 0, 4, 5], dtype=np.uint64), np.array([2, 0, 0, 0, 0]), np.array([0, 4, 4, 5, 5])), shape=(4,3))
obj.write_h5ad(noncanonical)
canonical = run('canonicalize-duplicates-zeros-order', noncanonical, plan)
assert json.loads((canonical / 'dataset.json').read_text())['matrix']['counts'] == [6,3,5]
for name, value in [('negative', -1.0), ('nan', float('nan')), ('infinity', float('inf')), ('fractional', 1.5), ('float-outside-exact-range', float(2**54))]:
    obj.layers['counts'] = sparse.csr_matrix(np.array([[value,0,0],[0,0,0],[0,0,0],[0,0,0]]))
    invalid = a.out / (name + '.h5ad')
    obj.write_h5ad(invalid)
    run('reject-' + name, invalid, plan, False)

run('reject-plan-version', source, dict(plan, schemaVersion=2), False)
run('reject-unknown-plan-field', source, dict(plan, sampleColum='sample'), False)
run('reject-normalized-X', source, dict(plan, matrixPath='X'), False)
run('reject-unknown-sample-column', source, dict(plan, sampleColumn='absent'), False)
run('reject-missing-design', source, dict(plan, sampleColumn='group'), False)
run('reject-invalid-layer-path', source, dict(plan, matrixPath='layers/../X'), False)

def malformed(name, mutate):
    path = a.out / (name + '.h5ad')
    shutil.copyfile(source, path)
    with h5py.File(path, 'r+') as f:
        mutate(f)
    run(name, path, plan, False)

malformed('reject-offsets', lambda f: f['layers/counts/indptr'].__setitem__(1, 1000))
malformed('reject-index', lambda f: f['layers/counts/indices'].__setitem__(0, 3))
malformed('reject-shape', lambda f: f['layers/counts'].attrs.__setitem__('shape', [4,4]))
malformed('reject-version', lambda f: f['layers/counts'].attrs.__setitem__('encoding-version', '9.0.0'))
malformed('reject-category', lambda f: f['obs/sample/codes'].__setitem__(0, -1))
def external(f):
    del f['layers/counts/data']
    f['layers/counts/data'] = h5py.ExternalLink(str(source.resolve()), '/layers/counts/data')
malformed('reject-external-link', external)

report = {'status':'passed-interoperability-not-biological-validation','checks':checks,'anndata':version('anndata'),'h5py':h5py.__version__, 'nativeHDF5':json.loads((a.out / 'csr' / 'runtime.json').read_text())['hdf5Version'], 'binarySHA256':hashlib.sha256(a.binary.read_bytes()).hexdigest()}
(a.out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
