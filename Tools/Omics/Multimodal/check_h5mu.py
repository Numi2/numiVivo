#!/usr/bin/env python3
"""Native H5MU import: complete measured CITE-seq plus explicit structural controls.

Requires outputs of check_citeseq.py. Source hashes and all command outcomes are
recorded. Structural spatial/ATAC checks do not qualify real spatial/ATAC assays.
"""
import argparse
import copy
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import time
from pathlib import Path
import h5py
import mudata
import numpy as np
from scipy import sparse


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def write(p, v):
    p.write_text(json.dumps(v, indent=2, allow_nan=False) + '\n')


p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--citeseq', type=Path, required=True)
p.add_argument('--plan', type=Path, default=Path(__file__).with_name('pbmc5k-h5mu-plan.json'))
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
commands = []


def run(label, args, success=True):
    start = time.monotonic()
    r = subprocess.run([str(a.binary), *map(str, args)], capture_output=True, text=True)
    (a.out / (label + '.log')).write_text(r.stdout + r.stderr)
    commands.append(dict(label=label, arguments=list(map(str, args)), expectedSuccess=success,
                         exitCode=r.returncode, seconds=time.monotonic() - start))
    write(a.out / 'commands.json', commands)
    assert (r.returncode == 0) == success, (label, r.stdout[-2000:], r.stderr[-2000:])


def ingest(label, source, plan, success=True):
    planpath = a.out / (label + '-plan.json')
    write(planpath, plan)
    dest = a.out / label
    run(label, ['multiassay-h5mu-import', source, '--plan', planpath, '--output', dest], success)
    if not success:
        assert not dest.exists(), label
        return
    assert sha(source) == sha(dest / 'original.h5')
    return dest


plan = json.loads(a.plan.read_text())
reference = a.citeseq / 'bundle'
assert sha(reference / 'original.h5') == '3b290ad9605b96974c9c16e5ae3427e5e5c496a66e55223df135b388b6d61417'
reference_hash = sha(reference / 'dataset.json')
# Full native CSR export, independent MuData rewrite, and mixed CSC/dense count
# layers must all produce exactly the previously qualified native dataset bytes.
sources = {'native-csr': reference / 'dataset.h5mu',
           'mudata-csr': a.citeseq / 'mudata-roundtrip.h5mu'}
mu = mudata.read_h5mu(sources['mudata-csr'])
mu.mod['rna'].layers['counts'] = mu.mod['rna'].X.tocsc()
mu.mod['protein'].layers['counts'] = mu.mod['protein'].X.toarray()
for mod in mu.mod.values():
    mod.X = sparse.csr_matrix(np.full((mod.n_obs, 1), .5)) @ sparse.csr_matrix(
        (np.ones(1), (np.zeros(1, dtype=int), np.zeros(1, dtype=int))), shape=(1, mod.n_vars))
sources['mixed-count-layers'] = a.out / 'mixed-count-layers.h5mu'
mu.write(sources['mixed-count-layers'])
del mu
for label, source in sources.items():
    selected = copy.deepcopy(plan)
    if label == 'mixed-count-layers':
        for assay in selected['assays']:
            assay['matrixPath'] = 'layers/counts'
    bundle = ingest(label, source, selected)
    assert sha(bundle / 'dataset.json') == reference_hash, label
    run(label + '-verify', ['multiassay-verify', bundle])
    with h5py.File(bundle / 'dataset.h5mu') as f, h5py.File(reference / 'dataset.h5mu') as ref:
        for key in ['rna', 'protein']:
            for col in ['data', 'indices', 'indptr']:
                np.testing.assert_array_equal(f[f'mod/{key}/X/{col}'][:], ref[f'mod/{key}/X/{col}'][:])
repeat = ingest('repeat', sources['mudata-csr'], plan)
for name in ['original.h5', 'plan.json', 'dataset.json', 'dataset.h5mu', 'receipt.json']:
    assert sha(repeat / name) == sha(a.out / 'mudata-csr' / name), name
run('overwrite', ['multiassay-h5mu-import', sources['mudata-csr'], '--plan', a.plan,
                  '--output', a.out / 'mudata-csr'], False)

# Small structural control: partial/reordered modality, exact UInt64, explicit
# spatial units, spot identity, whole-row missing positions and categorical obs.
fixture = json.loads((a.citeseq / 'fixture.json').read_text())
partial = mudata.read_h5mu(a.citeseq / 'partial.h5mu')
partial.obsm['spatial'] = np.array([[np.nan, np.nan], [np.nan, np.nan], [2.5, 4.]])
partial.write(a.out / 'partial-spatial.h5mu')
smallplan = copy.deepcopy(plan)
smallplan.update(id=fixture['id'], evidence='synthetic', sourceDescription=fixture['sourceDescription'],
                 samples=fixture['samples'], spatial=dict(path='obsm/spatial', frame=fixture['spatialFrames'][0]))
for mapped, native in zip(smallplan['assays'], fixture['assays']):
    for field in ['id', 'kind', 'featureNamespace', 'countUnit', 'genomeAssembly']:
        mapped[field] = native[field]
    # H5MU plan has one source description; identity/count comparison is separate.
    native['sourceDescription'] = fixture['sourceDescription']
small = ingest('partial-spatial', a.out / 'partial-spatial.h5mu', smallplan)
actual = json.loads((small / 'dataset.json').read_text())
# Swift omits absent optionals; normalize explicit JSON nulls for comparison.
def compact(v):
    if isinstance(v, dict):
        return {k: compact(x) for k, x in v.items() if x is not None}
    if isinstance(v, list):
        return [compact(x) for x in v]
    return v
assert actual == compact(fixture)
run('partial-verify', ['multiassay-verify', small])


def badfile(label, mutate):
    source = a.out / (label + '.h5mu')
    shutil.copy2(a.out / 'partial-spatial.h5mu', source)
    with h5py.File(source, 'r+') as f:
        mutate(f)
    ingest(label, source, smallplan, False)


def replace(f, path, values):
    attrs = dict(f[path].attrs)
    del f[path]
    d = f.create_dataset(path, data=values)
    d.attrs.update(attrs)


badfile('map-out-of-range', lambda f: f['obsmap/protein'].__setitem__(0, 9))
badfile('map-wrong-identity', lambda f: f['obsmap/protein'].__setitem__(slice(None), np.array([1, 0, 2]).reshape(f['obsmap/protein'].shape)))
badfile('map-omitted-local', lambda f: f['obsmap/protein'].__setitem__(0, 0))
badfile('varmap-unmapped', lambda f: f['varmap/protein'].__setitem__(slice(None), 0))
badfile('varmap-overlap', lambda f: f['varmap/protein'].__setitem__(slice(None), np.array([1, 0]).reshape(f['varmap/protein'].shape)))
badfile('axis-one', lambda f: f.attrs.__setitem__('axis', 1))
badfile('partial-nan', lambda f: f['obsm/spatial'].__setitem__((2, 0), np.nan))
badfile('infinite-position', lambda f: f['obsm/spatial'].__setitem__((2, 0), np.inf))
badfile('fractional-count', lambda f: replace(f, 'mod/protein/X/data', np.array([.5])))
badfile('negative-count', lambda f: replace(f, 'mod/protein/X/data', np.array([-1])))
# MuData writes repeated strings as categoricals. Change their category value.
def change_text(f, path, text):
    d = f[path]
    if isinstance(d, h5py.Group):
        d = d['categories'] if 'categories' in d else d['values']
    d[0] = text
for label, path, value in [
    ('sample-conflict', 'mod/protein/obs/sample', 'wrong-sample'),
    ('kind-conflict', 'mod/protein/obs/observation_kind', 'unknown-kind'),
    ('barcode-conflict', 'mod/protein/obs/barcode', 'wrong-barcode'),
    ('namespace-conflict', 'mod/protein/var/feature_namespace', 'wrong-namespace'),
    ('assay-kind-conflict', 'mod/protein/var/assay_kind', 'rna')]:
    badfile(label, lambda f, path=path, value=value: change_text(f, path, value))
for label, mutation in [
    ('omit-modality', lambda d: d.update(assays=d['assays'][:1])),
    ('unknown-option', lambda d: d.update(ignoreMissing=True)),
    ('missing-count-layer', lambda d: d['assays'][0].update(matrixPath='layers/missing')),
    ('raw-feature-axis', lambda d: d['assays'][0].update(matrixPath='raw/X')),
    ('unsafe-path', lambda d: d['assays'][0].update(matrixPath='layers/../X'))]:
    bad = copy.deepcopy(smallplan)
    mutation(bad)
    ingest(label, a.out / 'partial-spatial.h5mu', bad, False)
# Rehashed native data and native projection tampering must fail reconstruction.
for field in ['dataset', 'h5mu']:
    dest = a.out / ('tampered-' + field)
    shutil.copytree(small, dest)
    target = dest / ('dataset.json' if field == 'dataset' else 'dataset.h5mu')
    if field == 'dataset':
        changed = json.loads(target.read_text())
        changed['assays'][0]['matrix']['counts'][0] += 1
        target.write_text(json.dumps(changed, sort_keys=True, separators=(',', ':'), ensure_ascii=False))
    else:
        with h5py.File(target, 'r+') as f:
            f['mod/protein/X/data'][0] += np.uint64(1)
    receipt = json.loads((dest / 'receipt.json').read_text())
    assert field in receipt, receipt.keys()
    receipt[field] = {'bytes': list(bytes.fromhex(sha(target)))}
    write(dest / 'receipt.json', receipt)
    run('rehashed-' + field, ['multiassay-verify', dest], False)
write(a.out / 'checks.json', dict(status='passed', binarySHA256=sha(a.binary), checkerSHA256=sha(Path(__file__)),
    sources={k: dict(path=str(v), sha256=sha(v)) for k, v in sources.items()},
    referenceDatasetSHA256=reference_hash, completeRealCells=5247, rnaFeatures=33538, antibodyFeatures=32,
    fullDatasetExactAcrossLayouts=True, repeatedBundleBytesExact=True, partialMapSpatialUInt64Exact=True,
    commands=len(commands), expectedRejections=sum(not c['expectedSuccess'] for c in commands),
    packages={n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'mudata', 'anndata', 'h5py']},
    platform=platform.platform()))
