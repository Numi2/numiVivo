#!/usr/bin/env python3
"""Complete public paired RNA/ATAC counts: native 10x and independent MuData import.

This qualifies count/identity/interval interchange, not peak calling, biological
cell annotation, regulatory linkage or joint latent inference. No cell or feature
selection is allowed. ARC 2.0 peak matrix entries count cut sites, not fragments.
"""
import argparse
import copy
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import time
import anndata as ad
import h5py
import mudata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SHA = '5fbff5a4d85e0df345f6502e966ec787a8a4c429fd6b88a8772c43fd915cf3ff'
URL = 'https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5'
PAGE = 'https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0'

def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def write(p, obj):
    p.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--plan', type=Path, default=Path(__file__).with_name('pbmc3k-multiome-plan.json'))
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
assert sha(a.source) == SHA and a.source.stat().st_size == 38844318
a.out.mkdir(parents=True, exist_ok=False)
commands = []

def run(label, args, success=True):
    started = time.monotonic()
    r = subprocess.run(['/usr/bin/time', '-l', str(a.binary), *map(str, args)], capture_output=True, text=True)
    (a.out / (label + '.log')).write_text(r.stdout + r.stderr)
    rss = re.search(r'(\d+)\s+maximum resident set size', r.stderr)
    commands.append(dict(label=label, arguments=list(map(str,args)), exitCode=r.returncode,
                         expectedSuccess=success, seconds=time.monotonic()-started,
                         maximumResidentBytes=int(rss[1]) if rss else None))
    write(a.out / 'commands.json', commands)
    assert (r.returncode == 0) == success, (label, r.stdout[-1500:], r.stderr[-1500:])

plan = json.loads(a.plan.read_text())
bundle = a.out / 'bundle'
run('tenx-import', ['multiassay-10x-import', a.source, '--plan', a.plan, '--output', bundle])
run('tenx-verify', ['multiassay-verify', bundle])
run('tenx-repeat', ['multiassay-10x-import', a.source, '--plan', a.plan, '--output', a.out / 'repeat'])
for name in ['original.h5', 'plan.json', 'dataset.json', 'dataset.h5mu', 'receipt.json']:
    assert sha(bundle / name) == sha(a.out / 'repeat' / name), name
assert sha(bundle / 'original.h5') == SHA
# Independent direct source CSR construction and Scanpy reader must agree first.
with h5py.File(a.source) as f:
    m = f['matrix']
    original = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
    ids = m['features/id'].asstr()[:]
    names = m['features/name'].asstr()[:]
    kinds = m['features/feature_type'].asstr()[:]
    intervals = m['features/interval'].asstr()[:]
    barcodes = m['barcodes'].asstr()[:]
assert original.shape == (2711, 134920) and original.nnz == 24511186
scan = sc.read_10x_h5(a.source, gex_only=False)
assert (scan.X != original).nnz == 0
np.testing.assert_array_equal(scan.var['gene_ids'], ids)
np.testing.assert_array_equal(scan.obs_names, barcodes)
del scan
# Exact full native dataset, every count/axis and all genomic peak intervals.
native = json.loads((bundle / 'dataset.json').read_text())
assert len(native['observations']) == 2711 and native['samples'] == [plan['sample']]
assert [o['identity']['barcode'] for o in native['observations']] == list(barcodes)
mu = mudata.read_h5mu(bundle / 'dataset.h5mu')
reports = []
independent = {}
for key, kind, unit in [('rna', 'Gene Expression', 'umiCount'), ('atac', 'Peaks', 'cutSiteCount')]:
    mask = kinds == kind
    expected = original[:, mask]
    assay = next(v for v in native['assays'] if v['id'] == key)
    assert assay['observationIndices'] == list(range(2711)) and assay['countUnit'] == unit
    assert assay['genomeAssembly'] == 'GRCh38'
    assert [v['id'] for v in assay['features']] == list(ids[mask])
    assert [v['name'] for v in assay['features']] == list(names[mask])
    d = assay['matrix']
    count = sparse.csr_matrix((np.asarray(d['counts'], dtype=np.uint64), d['featureIndices'], d['rowOffsets']), shape=expected.shape)
    assert (count != expected).nnz == 0 and (mu.mod[key].X != expected).nnz == 0
    np.testing.assert_array_equal(mu.mod[key].var_names, ids[mask])
    np.testing.assert_array_equal(np.asarray(mu.obsmap[key]).ravel(), np.arange(1, 2712))
    np.testing.assert_array_equal(mu.mod[key].obs['barcode'], barcodes)
    if key == 'atac':
        for feature, source_interval in zip(assay['features'], intervals[mask]):
            contig, bounds = source_interval.split(':')
            start, end = map(int, bounds.split('-'))
            assert feature['interval'] == dict(contig=contig, start=start, end=end)
            assert feature['id'] == source_interval
    # Independent MuData is built from original source arrays, never native output.
    obs = pd.DataFrame(dict(sample=pd.Categorical([plan['sample']['id']] * 2711), barcode=barcodes,
                            observation_kind=pd.Categorical(['cell'] * 2711)), index=pd.Index(['observation-'+str(i) for i in range(2711)]))
    var = pd.DataFrame(dict(name=names[mask]), index=pd.Index(ids[mask]))
    independent[key] = ad.AnnData(X=expected.tocsc(), obs=obs, var=var)
    reports.append(dict(assay=key, unit=unit, cells=2711, features=expected.shape[1], nonzeros=expected.nnz,
                        totalCounts=int(expected.sum()), allCoordinatesExact=True, allFeatureIDsExact=True,
                        allPeakIntervalsExact=key == 'atac'))
    # Per-cell and per-feature totals independently check both axes.
    np.testing.assert_array_equal(np.asarray(count.sum(axis=0)), np.asarray(expected.sum(axis=0)))
    np.testing.assert_array_equal(np.asarray(count.sum(axis=1)), np.asarray(expected.sum(axis=1)))
    del count
assert reports[0]['totalCounts'] == 11786194 and reports[1]['totalCounts'] == 48245242
write(a.out / 'assays.json', reports)
del native, mu, original, assay, d
external = mudata.MuData(independent)
external.obs['sample'] = pd.Categorical([plan['sample']['id']] * 2711)
external.obs['barcode'] = barcodes
external.obs['observation_kind'] = pd.Categorical(['cell'] * 2711)
external.write(a.out / 'independent.h5mu')
del external, independent, expected
mapping = dict(schemaVersion=1, id=plan['id'], evidence=plan['evidence'], sourceDescription=plan['sourceDescription'],
               samples=[plan['sample']], sampleColumn='sample', barcodeColumn='barcode', defaultObservationKind='cell',
               observationKindColumn='observation_kind', assays=[])
for spec in plan['assays']:
    entry = {k: spec[k] for k in ['id', 'kind', 'featureNamespace', 'countUnit', 'genomeAssembly']}
    entry.update(sourceName=spec['id'], matrixPath='X', featureNameColumn='name')
    if spec['kind'] == 'chromatinAccessibility':
        entry['peakIDConvention'] = 'contig:start-end:zero-based-half-open'
    mapping['assays'].append(entry)
write(a.out / 'h5mu-plan.json', mapping)
run('independent-h5mu-import', ['multiassay-h5mu-import', a.out / 'independent.h5mu', '--plan', a.out / 'h5mu-plan.json', '--output', a.out / 'h5mu'])
run('independent-h5mu-verify', ['multiassay-verify', a.out / 'h5mu'])
assert sha(a.out / 'h5mu/dataset.json') == sha(bundle / 'dataset.json')
assert sha(a.out / 'h5mu/dataset.h5mu') == sha(bundle / 'dataset.h5mu')
# Invalid mapping and over-budget source reject before an artifact is published.
for label, change in [('missing-atac', lambda d: d.update(assays=d['assays'][:1])),
                      ('cutsite-rna', lambda d: d['assays'][0].update(countUnit='cutSiteCount')),
                      ('unknown-unit', lambda d: d['assays'][1].update(countUnit='accessibilityScore')),
                      ('assembly-mismatch', lambda d: d['assays'][1].update(genomeAssembly='GRCh37'))]:
    bad = copy.deepcopy(plan); change(bad); write(a.out / (label + '.json'), bad)
    run(label, ['multiassay-10x-import', a.source, '--plan', a.out / (label + '.json'), '--output', a.out / label], False)
    assert not (a.out / label).exists()
oversized = a.out / 'oversized.h5'
shutil.copy2(a.source, oversized)
with h5py.File(oversized, 'r+') as f:
    m = f['matrix']
    for col in ['data', 'indices']: del m[col]
    m.create_dataset('data', shape=(32_000_001,), dtype=np.uint64, chunks=(65536,))
    m.create_dataset('indices', shape=(32_000_001,), dtype=np.uint32, chunks=(65536,))
run('entry-budget', ['multiassay-10x-import', oversized, '--plan', a.plan, '--output', a.out / 'oversized'], False)
assert not (a.out / 'oversized').exists()
write(a.out / 'checks.json', dict(status='passed', sourceURL=URL, sourcePage=PAGE, sourceSHA256=SHA,
    binarySHA256=sha(a.binary), checkerSHA256=sha(Path(__file__)), planSHA256=sha(a.plan),
    independentH5MUSHA256=sha(a.out/'independent.h5mu'), independentDatasetBytesExact=True,
    repeatedBundleBytesExact=True, allSourceCellsAndFeaturesRetained=True, assays=reports,
    commands=len(commands), expectedRejections=sum(not c['expectedSuccess'] for c in commands),
    measuredPeakResidentBytes=max(c['maximumResidentBytes'] or 0 for c in commands),
    scope='Complete-source count, pairing and interval interchange; not biological ATAC/joint-model qualification',
    platform=platform.platform(), packages={n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'anndata', 'mudata', 'scanpy', 'h5py']}))
