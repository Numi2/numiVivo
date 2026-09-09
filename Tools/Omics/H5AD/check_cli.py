#!/usr/bin/env python3
"""Check the full product H5AD commands and their existing count receipt route."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import anndata
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--plan', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
checks = []

def run(name, *args, success=True):
    result = subprocess.run([str(a.binary), *map(str, args)], capture_output=True, text=True)
    (a.out / (name + '.log')).write_text(result.stdout + result.stderr)
    assert (result.returncode == 0) == success, (name, result.stderr)
    checks.append(name)
    return json.loads(result.stdout) if success and result.stdout.strip() else None

imported = a.out / 'imported'
run('import', 'singlecell-h5ad-import', a.source, '--plan', a.plan, '--output', imported)
assert (imported / 'original.h5ad').read_bytes() == a.source.read_bytes()
receipt = json.loads((imported / 'h5ad-import.json').read_text())
# Verify the existing fingerprint byte-array encoding.
for key, name in [('source','original.h5ad'),('dataset','dataset.json'),('mapping','h5ad-mapping.json')]:
    assert bytes(receipt[key]['bytes']) == hashlib.sha256((imported / name).read_bytes()).digest()
run('refuse-import-overwrite', 'singlecell-h5ad-import', a.source, '--plan', a.plan, '--output', imported, success=False)
exported = a.out / 'exported.h5ad'
run('export', 'singlecell-h5ad-write', imported / 'dataset.json', '--output', exported)
run('refuse-export-overwrite', 'singlecell-h5ad-write', imported / 'dataset.json', '--output', exported, success=False)
actual = anndata.read_h5ad(exported)
expected = anndata.read_h5ad(a.source).layers['counts']
assert np.array_equal(actual.X.toarray(), expected.toarray())
counts_receipt = a.out / 'counts-receipt.json'
run('count-campaign', 'singlecell-run', imported / 'manifest.json', '--store', a.out / 'store', '--output', counts_receipt)
run('count-replay', 'singlecell-verify', counts_receipt, '--store', a.out / 'store')
report = dict(status='passed-full-product-H5AD-CLI', checks=checks, binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(), hdf5Version=receipt['hdf5Version'])
(a.out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
