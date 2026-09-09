#!/usr/bin/env python3
"""Qualify source-preserving native edits against the actual AnnData reader."""
import argparse
import copy
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--full-product', action='store_true')
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
checks = []

def fingerprint(path):
    return {'bytes': list(hashlib.sha256(path.read_bytes()).digest())}

def val(kind, **fields):
    return {kind: fields}

def edit(path, value, mode='add'):
    return dict(path=path, mode=mode, value=value)

def run(name, source, plan, success=True, destination=None):
    before = source.read_bytes()
    planpath = a.out / (name + '.json')
    planpath.write_text(json.dumps(plan, ensure_ascii=False))
    destination = destination or a.out / (name + '.h5ad')
    old_output = destination.read_bytes() if destination.exists() else None
    command = ([str(a.binary), 'singlecell-h5ad-annotate', str(source), '--plan', str(planpath), '--output', str(destination)] if a.full_product
               else [str(a.binary), 'annotate', str(source), str(planpath), str(destination)])
    result = subprocess.run(command, capture_output=True, text=True)
    (a.out / (name + '.log')).write_text(result.stdout + result.stderr)
    assert result.returncode == (0 if success else 65), (name, result.returncode, result.stderr)
    assert source.read_bytes() == before, 'source was modified'
    if success:
        receipt = json.loads(result.stdout)
        assert receipt['source'] == fingerprint(source)
        assert receipt['output'] == fingerprint(destination)
        checks.append(name)
        return destination, receipt
    assert (destination.read_bytes() if destination.exists() else None) == old_output, 'failed operation published or modified output'
    checks.append(name)

# Source has measured-looking structure but is explicitly a software fixture.
obj = ad.AnnData(X=sparse.csr_matrix(np.array([[1,0,0],[0,2,0],[0,0,3],[4,0,0]], dtype=np.uint64)),
    obs=pd.DataFrame({'sample': pd.Categorical(['s1','s1','s2','s2'], ordered=True), 'old_nullable': pd.array([1,None,2,3], dtype='Int64')}, index=['c1','c2','c3','c4']),
    var=pd.DataFrame({'name': ['α','b','c']}, index=['g1','g2','g3']))
obj.layers['old_counts'] = obj.X.copy()
obj.obsm['old_embedding'] = np.arange(8, dtype=np.float32).reshape(4,2)
obj.varm['old_loadings'] = np.arange(6, dtype=np.float32).reshape(3,2)
obj.obsp['old_graph'] = sparse.eye(4, format='csc', dtype=np.float32)
obj.uns['old_nested'] = {'labels': np.array(['α','β'], dtype=object), 'precision': np.float32(0.25)}
obj.raw = ad.AnnData(X=sparse.csr_matrix(np.arange(20, dtype=np.int64).reshape(4,5)), obs=obj.obs.copy(), var=pd.DataFrame(index=['r1','r2','r3','r4','r5']))
source = a.out / 'source.h5ad'
obj.write_h5ad(source)
with h5py.File(source, 'r+') as f:
    # Editing obs must not change an existing hard-linked backup in uns.
    f['uns']['obs_backup'] = f['obs']
    f['uns']['embedding_backup'] = f['obsm/old_embedding']
plan = dict(schemaVersion=1, source=fingerprint(source), provenance='Synthetic interoperability fixture; annotations are not biological claims', edits=[
    edit('obs/score', val('float64', shape=[4], values=[-1,0,0.5,1])),
    edit('obs/cluster', val('categorical', codes=[1,-1,0,1], categories=['B','T','unused'], ordered=True)),
    edit('obs/nullable_float', val('nullableFloat64', values=[0.5,None,-0.5,0])),
    edit('obs/nullable_integer', val('nullableInt64', values=[9007199254740993,None,-4,0])),
    edit('obs/nullable_boolean', val('nullableBoolean', values=[True,None,False,True])),
    edit('obs/nullable_string', val('nullableString', values=['λ',None,'','a'])),
    edit('obs/int64', val('int64', shape=[4], values=[-(2**63),2**63-1,-1,0])),
    edit('obs/uint64', val('uint64', shape=[4], values=[2**64-1,0,9007199254740993,4])),
    edit('var/highly_variable', val('boolean', shape=[3], values=[True,False,True])),
    edit('obsm/native_pca', val('float64', shape=[4,2], values=list(range(8)))),
    edit('varm/native_loadings', val('float64', shape=[3,2], values=list(range(6)))),
    edit('obsp/native_distances', val('csrFloat64', shape=[4,4], indptr=[0,1,2,3,4], indices=[1,2,3,0], values=[0.25,0.5,0.75,1])),
    edit('varp/native_graph', val('csrUInt64', shape=[3,3], indptr=[0,1,2,3], indices=[0,1,2], values=[1,1,1])),
    edit('layers/native_values', val('csrFloat64', shape=[4,3], indptr=[0,1,2,3,4], indices=[0,1,2,0], values=[0.25,0.5,0.75,1])),
    edit('uns/native_results', val('dictionary', values={
        'description': val('string', shape=[], values=['Native test result']),
        'labels': val('string', shape=[2,2], values=['a','β','c','d']),
        'parameter': val('float64', shape=[], values=[0.5]),
        'enabled': val('boolean', shape=[], values=[True]),
        'count': val('uint64', shape=[], values=[2**64-1])
    }))
])
output, receipt = run('typed-annotations', source, plan)
actual = ad.read_h5ad(output)
assert actual.obs['uint64'].iloc[0] == 2**64-1
assert actual.obs['int64'].iloc[0] == -(2**63)
assert actual.obs['nullable_integer'].iloc[0] == 9007199254740993
assert pd.isna(actual.obs['nullable_integer'].iloc[1])
assert np.isnan(actual.obs['nullable_float'].iloc[1])
assert actual.obs['nullable_float'].iloc[2] == -0.5
assert actual.obs['nullable_string'].iloc[2] == ''
assert pd.isna(actual.obs['cluster'].iloc[1])
assert list(actual.obs['cluster'].cat.categories) == ['B','T','unused'] and actual.obs['cluster'].cat.ordered
assert actual.var['highly_variable'].tolist() == [True,False,True]
assert np.array_equal(actual.obsm['native_pca'], np.arange(8).reshape(4,2))
assert actual.obsp['native_distances'][0,1] == 0.25
assert actual.layers['native_values'][3,0] == 1
assert actual.uns['native_results']['labels'].shape == (2,2)
assert actual.uns['native_results']['count'] == 2**64-1
assert 'score' not in actual.uns['obs_backup'].columns
assert actual.raw.shape == (4,5)
# Check EVERY existing stored dataset and its metadata, not just selected fields.
with h5py.File(source) as old, h5py.File(output) as new:
    def compare(name, node):
        assert name in new
        for key, value in node.attrs.items():
            if key == 'column-order' and name in ['obs','var']:
                assert list(new[name].attrs[key])[:len(value)] == list(value)
            else:
                np.testing.assert_array_equal(new[name].attrs[key], value)
        if isinstance(node, h5py.Dataset):
            assert node.dtype == new[name].dtype and node.shape == new[name].shape
            np.testing.assert_array_equal(node[()], new[name][()])
    old.visititems(compare)
checks.append('all-untouched-datasets-dtypes-attributes-preserved')
history = actual.uns['numivivo_edits']
assert len(history) == 1
record = next(iter(history.values()))
assert record['source_sha256'] == bytes(receipt['source']['bytes']).hex()
assert record['plan_sha256'] == bytes(receipt['plan']['bytes']).hex()
assert json.loads(record['plan_json']) == plan
checks.append('embedded-plan-and-source-provenance')
second_plan = dict(schemaVersion=1, source=fingerprint(output), provenance='Explicit second test edit', edits=[
    edit('obs/score', val('float64', shape=[4], values=[4,3,2,1]), mode='replace'),
    edit('obsm/old_embedding', val('float64', shape=[4,2], values=[0]*8), mode='replace')])
second, _ = run('explicit-replacement-retains-history', output, second_plan)
second_obj = ad.read_h5ad(second)
assert second_obj.obs['score'].tolist() == [4,3,2,1]
assert np.array_equal(second_obj.uns['embedding_backup'], obj.obsm['old_embedding'])
assert len(second_obj.uns['numivivo_edits']) == 2
run('refuse-existing-output', source, plan, False, output)

def reject(name, **changes):
    invalid = copy.deepcopy(plan)
    invalid.update(changes)
    run(name, source, invalid, False)
reject('reject-source-mismatch', source={'bytes':[0]*32})
reject('reject-plan-version', schemaVersion=2)
reject('reject-unknown-plan-field', unknown=True)
reject('reject-missing-provenance', provenance=' ')
reject('reject-column-length', edits=[edit('obs/bad', val('int64', shape=[3], values=[1,2,3]))])
reject('reject-index-replacement', edits=[edit('obs/_index', val('string', shape=[4], values=['x']*4), 'replace')])
reject('reject-duplicate-edits', edits=[plan['edits'][0], plan['edits'][0]])
reject('reject-shape-overflow', edits=[edit('uns/bad', val('float64', shape=[2**62,2**62], values=[]))])
reject('reject-embedding-axis', edits=[edit('obsm/bad', val('float64', shape=[2,2], values=[0]*4))])
reject('reject-dense-layer', edits=[edit('layers/bad', val('float64', shape=[4,3], values=[0]*12))])
reject('reject-invalid-CSR', edits=[edit('obsp/bad', val('csrFloat64', shape=[4,4], indptr=[0,1,2,3,5], indices=[0,1,2,3], values=[1]*4))])
reject('reject-category-range', edits=[edit('obs/bad', val('categorical', codes=[0,0,0,2], categories=['only'], ordered=False))])
reject('reject-reserved-history', edits=[edit('uns/numivivo_edits', val('dictionary', values={}))])
reject('reject-path-traversal', edits=[edit('obs/../uns', val('string', shape=[], values=['bad']))])
reject('failed-late-edit-publishes-nothing', edits=[plan['edits'][0], edit('uns/nonexistent', val('int64', shape=[], values=[1]), mode='replace')])
reject('reject-add-existing-field', edits=[edit('obs/sample', val('string', shape=[4], values=['x']*4))])
# Storage dependencies must not silently change meaning when moved to a new file.
for kind in ['external-link','soft-link','virtual-dataset','object-reference','attribute-reference']:
    variant = a.out / (kind + '-source.h5ad')
    variant.write_bytes(source.read_bytes())
    with h5py.File(variant, 'r+') as f:
        if kind == 'external-link':
            f['uns/external'] = h5py.ExternalLink(str(source.resolve()), '/X')
        elif kind == 'soft-link':
            f['uns/soft'] = h5py.SoftLink('/X')
        elif kind == 'virtual-dataset':
            layout = h5py.VirtualLayout(shape=(4,2), dtype=np.float32)
            layout[:] = h5py.VirtualSource(str(source.resolve()), '/obsm/old_embedding', shape=(4,2))
            f['uns'].create_virtual_dataset('virtual', layout)
        elif kind == 'object-reference':
            f['uns'].create_dataset('ref', data=np.array([f['X'].ref], dtype=h5py.ref_dtype))
        else:
            f['obs'].attrs['reference'] = f['X'].ref
    changed = copy.deepcopy(plan)
    changed['source'] = fingerprint(variant)
    run('reject-' + kind, variant, changed, False)
report = dict(status='passed-native-AnnData-annotation-interoperability-not-biological-validation', checks=checks,
    anndata=version('anndata'), h5py=h5py.__version__, fullProduct=a.full_product,
    binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(), nativeHDF5=receipt['hdf5Version'])
(a.out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
