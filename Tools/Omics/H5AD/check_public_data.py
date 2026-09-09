#!/usr/bin/env python3
"""Full PBMC3K annotation preservation, not a DE/integration benchmark.

Download URL and checksum are pinned. No cells/features are selected or dropped.
The upstream legacy H5AD is retained and re-encoded by the current reference
AnnData writer before entering the native versioned H5AD path.
"""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import urllib.request
import anndata as ad
import h5py
import numpy as np
import pandas as pd

URL = 'https://falexwolf.de/data/pbmc3k_raw.h5ad'
SHA256 = '89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1'
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--download', type=Path, help='Use an existing downloaded file, still verified against the pinned checksum')
p.add_argument('--full-product', action='store_true')
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
raw = a.out / 'pbmc3k_raw.h5ad'
if a.download:
    raw.write_bytes(a.download.read_bytes())
else:
    with urllib.request.urlopen(URL, timeout=120) as response:
        data = response.read(8 * 1024 * 1024 + 1)
    assert len(data) <= 8 * 1024 * 1024
    raw.write_bytes(data)
assert hashlib.sha256(raw.read_bytes()).hexdigest() == SHA256
original = ad.read_h5ad(raw)
assert original.shape == (2700,32738) and original.X.nnz == 2286884
modern = a.out / 'reference-current.h5ad'
original.write_h5ad(modern, compression='gzip')
reference = ad.read_h5ad(modern)
np.testing.assert_array_equal(reference.X.data, original.X.data)
np.testing.assert_array_equal(reference.X.indices, original.X.indices)
np.testing.assert_array_equal(reference.X.indptr, original.X.indptr)
assert reference.X.dtype == original.X.dtype
np.testing.assert_array_equal(reference.obs_names, original.obs_names)
np.testing.assert_array_equal(reference.var_names, original.var_names)
plan = dict(schemaVersion=1, source={'bytes':list(hashlib.sha256(modern.read_bytes()).digest())},
    provenance='Public PBMC3K file preservation check; original download SHA256=' + SHA256 + '; reference AnnData re-encoding=' + version('anndata') + '; no cell selection, biological labeling or statistical fit',
    edits=[dict(path='obs/numivivo_source_library', mode='add', value={'string':{'shape':[2700],'values':['pbmc3k-public-library']*2700}}),
           dict(path='uns/numivivo_interop', mode='add', value={'dictionary':{'values':{
               'download_url':{'string':{'shape':[],'values':[URL]}},
               'download_sha256':{'string':{'shape':[],'values':[SHA256]}},
               'reference_anndata':{'string':{'shape':[],'values':[version('anndata')]}},
               'biological_validation':{'boolean':{'shape':[],'values':[False]}}
           }}})])
planpath = a.out / 'annotations.json'
planpath.write_text(json.dumps(plan))
output = a.out / 'native-annotated.h5ad'
command = ([str(a.binary), 'singlecell-h5ad-annotate', str(modern), '--plan', str(planpath), '--output', str(output)] if a.full_product
    else [str(a.binary), 'annotate', str(modern), str(planpath), str(output)])
result = subprocess.run(command, capture_output=True, text=True)
(a.out / 'native.log').write_text(result.stdout + result.stderr)
assert result.returncode == 0, result.stderr
receipt = json.loads(result.stdout)
assert bytes(receipt['output']['bytes']) == hashlib.sha256(output.read_bytes()).digest()
restored = ad.read_h5ad(output)
assert restored.shape == original.shape and restored.X.dtype == original.X.dtype
np.testing.assert_array_equal(restored.X.data, original.X.data)
np.testing.assert_array_equal(restored.X.indices, original.X.indices)
np.testing.assert_array_equal(restored.X.indptr, original.X.indptr)
np.testing.assert_array_equal(restored.obs_names, original.obs_names)
pd.testing.assert_frame_equal(restored.var, reference.var)
pd.testing.assert_frame_equal(restored.obs.drop(columns=['numivivo_source_library']), reference.obs)
assert restored.obs['numivivo_source_library'].eq('pbmc3k-public-library').all()
assert not restored.uns['numivivo_interop']['biological_validation']
# Also compare every original stored dataset/dtype/attribute after re-encoding.
with h5py.File(modern) as source, h5py.File(output) as target:
    def verify(name, node):
        for key, value in node.attrs.items():
            if name == 'obs' and key == 'column-order':
                assert list(target[name].attrs[key])[:len(value)] == list(value)
            else:
                np.testing.assert_array_equal(target[name].attrs[key], value)
        if isinstance(node, h5py.Dataset):
            assert target[name].dtype == node.dtype and target[name].shape == node.shape
            np.testing.assert_array_equal(target[name][()], node[()])
    source.visititems(verify)
report = dict(status='passed-full-public-PBMC3K-annotation-preservation', sourceURL=URL, sourceSHA256=SHA256,
    cells=2700, features=32738, nonzeros=2286884, dtype=str(original.X.dtype), cellOrFeatureSubsetting=False,
    referenceReencodingRequired=True, directLegacyH5ADSupport=False, countProjectionQualified=False,
    countProjectionLimitation='2286884 nonzeros exceeds current default count-model limit of 2000000; this run qualifies source-preserving annotations only',
    biologicalBenchmark=False, annotationReceipt=receipt, anndata=version('anndata'), h5py=h5py.__version__,
    binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest())
(a.out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k:v for k,v in report.items() if k != 'annotationReceipt'}, indent=2))
