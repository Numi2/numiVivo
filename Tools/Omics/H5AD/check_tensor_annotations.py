#!/usr/bin/env python3
"""Native tensor writes and projection checked against AnnData; no biology claims."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--source', type=Path)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
fingerprint = lambda p: {'bytes': list(bytes.fromhex(sha(p)))}
obj = ad.read_h5ad(a.source) if a.source else ad.AnnData(
    X=sparse.eye(4, 3, format='csr'),
    obs=pd.DataFrame(index=['c0', 'c1', 'c2', 'c3']),
    var=pd.DataFrame(index=['g0', 'g1', 'g2']))
source = a.out / 'source.h5ad'
obj.write_h5ad(source)
original_hash = sha(source)
cells, features = obj.shape
arrays = {
    'obsm/tensor': np.arange(cells * 4, dtype=np.float64).reshape(cells, 2, 2),
    'varm/tensor': np.arange(features * 2, dtype=np.int64).reshape(features, 1, 2),
    'obsm/rank8': np.arange(cells, dtype=np.uint64).reshape(cells, 1, 1, 1, 1, 1, 1, 1),
    'obsm/empty': np.zeros((cells, 2, 0), dtype=np.float64),
    'uns/boolean_tensor': np.ones((2, 1, 2), dtype=bool),
    'uns/string_tensor': np.array(['a', 'β', '', 'd']).reshape(2, 1, 2),
}
def edit(path, x):
    kind = {'f': 'float64', 'i': 'int64', 'u': 'uint64', 'b': 'boolean', 'U': 'string'}[x.dtype.kind]
    return {'path': path, 'mode': 'add', 'value': {kind: {'shape': list(x.shape), 'values': x.ravel().tolist()}}}
plan = {'schemaVersion': 1, 'source': fingerprint(source),
        'provenance': 'Tensor interoperability control; derived values are not biological predictions',
        'edits': [edit(path, x) for path, x in arrays.items()]}
checks = []
def run(tag, args, success=True):
    result = subprocess.run([str(a.binary), *args], capture_output=True, text=True)
    (a.out / (tag + '.log')).write_text(result.stdout + result.stderr)
    assert result.returncode == (0 if success else 65), (tag, result.returncode, result.stderr)
    assert sha(source) == original_hash
    checks.append(tag)
    return result
planpath = a.out / 'plan.json'
planpath.write_text(json.dumps(plan))
output = a.out / 'annotated.h5ad'
run('annotate', ['singlecell-h5ad-annotate', str(source), '--plan', str(planpath), '--output', str(output)])
got = ad.read_h5ad(output)
for path, expected in arrays.items():
    slot, key = path.split('/')
    np.testing.assert_array_equal(getattr(got, slot)[key], expected)
    if expected.dtype.kind != 'U':
        assert getattr(got, slot)[key].dtype == expected.dtype
assert (got.X != obj.X).nnz == 0
pd.testing.assert_frame_equal(got.obs, obj.obs)
pd.testing.assert_frame_equal(got.var, obj.var)
rows, columns = [cells - 1, 0, cells - 1], [features - 1, 0]
projection = {'schemaVersion': 1, 'source': fingerprint(output),
              'provenance': 'Repeated tensor-axis interoperability control',
              'observationIndices': rows, 'featureIndices': columns}
pp = a.out / 'projection.json'
pp.write_text(json.dumps(projection))
bundle = a.out / 'projection'
run('project', ['singlecell-h5ad-project', str(output), '--plan', str(pp), '--output', str(bundle)])
projected = ad.read_h5ad(bundle / 'projected.h5ad')
for path, expected in arrays.items():
    slot, key = path.split('/')
    if slot == 'obsm': expected = expected[rows]
    if slot == 'varm': expected = expected[columns]
    np.testing.assert_array_equal(getattr(projected, slot)[key], expected)
run('replay', ['singlecell-h5ad-project-verify', str(bundle)])
for tag, path, shape, values in [
    ('rank9', 'uns/bad', [1] * 9, [0]),
    ('wrong-axis', 'obsm/bad', [cells + 1, 1, 1], [0] * (cells + 1)),
    ('components', 'obsm/bad', [cells, 257, 0], []),
    ('shape-mismatch', 'uns/bad', [2, 2, 2], [0]),
    ('overflow', 'uns/bad', [2**62, 2, 2], []),
]:
    bad = dict(plan, edits=[{'path': path, 'mode': 'add', 'value': {'float64': {'shape': shape, 'values': values}}}])
    bp, dest = a.out / (tag + '.json'), a.out / (tag + '.h5ad')
    bp.write_text(json.dumps(bad))
    run(tag, ['singlecell-h5ad-annotate', str(source), '--plan', str(bp), '--output', str(dest)], False)
    assert not dest.exists()
report = {'status': 'passed-software-interoperability', 'checks': checks,
          'shape': list(obj.shape), 'anndataVersion': ad.__version__,
          'binarySHA256': sha(a.binary), 'checkerSHA256': sha(Path(__file__)),
          'sourceSHA256': original_hash, 'outputSHA256': sha(output),
          'upstreamSourceSHA256': sha(a.source) if a.source else None}
(a.out / 'report.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report))
