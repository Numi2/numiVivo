#!/usr/bin/env python3
"""Verify the complete prepared source through AnnData's backed public reader."""
import argparse
import json
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np

from acquire import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    source = root / 'hirisa.h5ad'
    prepared = json.loads((root / 'prepared-receipt.json').read_text())
    receipt = json.loads((root / 'source-receipt.json').read_text())
    summary = json.loads((root / 'source-audit/summary.json').read_text())
    assert prepared['complete'] and summary['complete']
    assert digest(source) == prepared['sourceSHA256']
    data = ad.read_h5ad(source, backed='r')
    try:
        assert data.isbacked and data.shape == (summary['cells'], summary['features'])
        assert data.obs_names.is_unique and data.var_names.is_unique
        assert np.array_equal(data.obs_names, data.obs['cell_uuid'].astype(str))
        assert set(data.obs.geo_accession.astype(str)) == {s['accession'] for s in receipt['files']}
        offset = 0
        nnz = total = 0
        for file in receipt['files']:
            accession = file['accession']
            audit = json.loads((root / 'source-audit' / (accession + '.json')).read_text())
            n = audit['cells']
            assert data.obs.geo_accession.iloc[offset:offset + n].astype(str).eq(accession).all()
            reference_path = root / 'source-audit' / audit['referencePath']
            assert digest(reference_path) == audit['referenceSHA256']
            with np.load(reference_path) as reference:
                gene_counts = np.zeros(data.n_vars, dtype=np.uint64)
                for start in range(0, n, 1024):
                    stop = min(start + 1024, n)
                    block = data.X[offset + start:offset + stop].astype(np.uint64)
                    block.sum_duplicates()
                    block.eliminate_zeros()
                    row_counts = np.asarray(block.sum(axis=1)).ravel()
                    assert np.array_equal(row_counts, reference['totalCounts'][start:stop])
                    assert np.array_equal(np.diff(block.indptr), reference['detectedFeatures'][start:stop])
                    gene_counts += np.asarray(block.sum(axis=0)).ravel()
                    nnz += block.nnz
                    total += int(row_counts.sum())
                assert np.array_equal(gene_counts, reference['geneCounts'])
            offset += n
            print(json.dumps({'accession': accession, 'verifiedCells': offset}), flush=True)
        assert offset == summary['cells'] and nnz == summary['canonicalNonzeros'] and total == summary['totalUMIs']
        result = {'status': 'passed', 'sourceSHA256': prepared['sourceSHA256'],
                  'verifierSHA256': digest(Path(__file__)), 'anndataVersion': version('anndata'),
                  'cells': offset, 'features': data.n_vars, 'canonicalNonzeros': nnz, 'UMIs': total,
                  'sourceLibraryAggregatesVerified': len(receipt['files']),
                  'allCellQCVerified': True, 'backedSparseMatrix': True,
                  'observationColumns': list(data.obs.columns),
                  'nativeExecutionQualified': False, 'biologicalModelValidation': False}
        (root / 'anndata-verification.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
        print(json.dumps(result), flush=True)
    finally:
        data.file.close()

if __name__ == '__main__':
    main()
