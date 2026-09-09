#!/usr/bin/env python3
"""Full original Baron rows -> sparse AnnData -> native streamed QC/pseudobulk."""
import argparse
import csv
import gzip
import hashlib
from importlib.metadata import version
import io
import json
from pathlib import Path
import tarfile

import anndata as ad
import numpy as np
import scanpy as sc

from prepare_baron import ARCHIVE_SHA, digest

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--archive', type=Path, required=True)
p.add_argument('--prepared', type=Path, required=True)
p.add_argument('--bundle', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
assert digest(a.archive) == ARCHIVE_SHA
obj = ad.read_h5ad(a.prepared/'prepared.h5ad')
report = json.loads((a.bundle/'report.json').read_text())
receipt = json.loads((a.bundle/'receipt.json').read_text())
assert obj.shape == (8569, 20125) and obj.X.nnz == 16171764
assert report['canonicalNonzeros'] == obj.X.nnz and report['contrasts'] == []
assert [f['id'] for f in report['metadata']['features']] == obj.var_names.tolist()
assert [c['barcode'] for c in report['metadata']['cells']] == obj.obs_names.tolist()
assert [c['group'] for c in report['metadata']['cells']] == obj.obs.cell_type.tolist()
assert [c['sampleID'] for c in report['metadata']['cells']] == obj.obs.native_sample.tolist()

# Verify every stored row against the original author CSV, independent of the
# native scanner. A single row is dense; no cells-by-genes array is allocated.
row_number = 0
with tarfile.open(a.archive) as archive:
    for donor in range(1, 5):
        name = f'GSM{2230756+donor}_human{donor}_umifm_counts.csv.gz'
        with gzip.GzipFile(fileobj=archive.extractfile(name)) as zipped:
            rows = csv.reader(io.TextIOWrapper(zipped))
            header = next(rows)
            assert header[3:] == obj.var_names.tolist()
            for row in rows:
                original = np.asarray(row[3:], dtype=np.int64)
                assert row[0] == obj.obs_names[row_number]
                assert row[1] == obj.obs.source_barcode.iloc[row_number]
                assert row[2] == obj.obs.cell_type.iloc[row_number]
                lo, hi = obj.X.indptr[row_number:row_number+2]
                indices = np.flatnonzero(original)
                assert np.array_equal(indices, obj.X.indices[lo:hi])
                assert np.array_equal(original[indices], obj.X.data[lo:hi])
                assert int(original.sum()) == report['quality'][row_number]['totalCounts']
                assert len(indices) == report['quality'][row_number]['detectedFeatures']
                row_number += 1
assert row_number == len(obj)
qc, _ = sc.pp.calculate_qc_metrics(obj, percent_top=None, log1p=False, inplace=False)
assert np.array_equal(qc.total_counts.to_numpy(), [q['totalCounts'] for q in report['quality']])
assert np.array_equal(qc.n_genes_by_counts.to_numpy(), [q['detectedFeatures'] for q in report['quality']])
assert all(q.get('mitochondrialFraction') is None for q in report['quality'])

bulk = report['pseudobulk']
matrix = bulk['matrix']
membership = []
assert len(bulk['groups']) == 56 and bulk['featureIDs'] == obj.var_names.tolist()
for i, group in enumerate(bulk['groups']):
    assert len(group['sampleIDs']) == 1
    donor = group['sampleIDs'][0]
    assert group['donorID'] == donor
    assert group['condition'] == ('T2D' if donor == 'human4' else 'nonT2D')
    mask = (obj.obs.native_sample == donor) & (obj.obs.cell_type == group['cellGroup'])
    selected = np.flatnonzero(mask)
    assert group['sourceCellIndices'] == selected.tolist()
    membership.extend(selected.tolist())
    expected = np.asarray(obj.X[selected].sum(axis=0)).ravel()
    lo, hi = matrix['rowOffsets'][i:i+2]
    nonzero = np.flatnonzero(expected)
    assert nonzero.tolist() == matrix['featureIndices'][lo:hi]
    assert expected[nonzero].tolist() == matrix['counts'][lo:hi]
assert sorted(membership) == list(range(len(obj)))
source_sha = digest(a.bundle/'original.h5ad')
assert source_sha == digest(a.prepared/'prepared.h5ad')
assert list(bytes.fromhex(source_sha)) == receipt['source']['bytes']
result = dict(status='passed-full-source-QC-and-pseudobulk', cells=len(obj), features=obj.n_vars,
    nonzeros=obj.X.nnz, totalUMIs=int(obj.X.sum()), donors=4, cellTypes=14, pseudobulks=56,
    archiveSHA256=ARCHIVE_SHA, preparedSHA256=source_sha, reportSHA256=digest(a.bundle/'report.json'),
    receipt=receipt, checks=['Every original author CSV row exactly matches sparse H5AD and native total/detected counts',
        'Scanpy total/detected counts exactly match all native cells', 'All 56 donor/cell-type pseudobulks equal sparse source sums',
        'Every source cell occurs exactly once in pseudobulk membership', 'All source feature identities/order and cell labels retained',
        'T2D donor status preserved; missing mitochondrial annotation remains missing'],
    versions={name: version(name) for name in ['anndata', 'numpy', 'scipy', 'scanpy']},
    qualification='Full-source count/QC/pseudobulk evidence only. Not donor integration, replicated disease DE, cell-type ground truth, or million-cell/out-of-core PCA qualification.')
a.out.parent.mkdir(parents=True, exist_ok=True)
a.out.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({key: result[key] for key in ['status', 'cells', 'features', 'nonzeros', 'totalUMIs', 'pseudobulks']}))
