#!/usr/bin/env python3
"""Reconstruct mitochondrial QC from the complete backed sparse source."""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
from acquire import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    a = p.parse_args();root = a.root
    prepared = json.loads((root / 'prepared-receipt.json').read_text())
    source = json.loads((root / 'source-receipt.json').read_text())
    mapping = json.loads((root / 'stream-plan.json').read_text())['mapping']
    assert prepared['complete'] and digest(root / 'hirisa.h5ad') == prepared['sourceSHA256']
    output = root / 'native-qc-reference';output.mkdir(exist_ok=False)
    data = ad.read_h5ad(root / 'hirisa.h5ad', backed='r')
    records = [];offset = 0
    try:
        mito = np.flatnonzero(data.var_names.isin(mapping['mitochondrialFeatureIDs']))
        assert len(mito) == len(set(mapping['mitochondrialFeatureIDs']))
        for file in source['files']:
            accession = file['accession']
            audit = json.loads((root / 'source-audit' / (accession + '.json')).read_text())
            reference = root / 'source-audit' / audit['referencePath']
            assert digest(reference) == audit['referenceSHA256']
            n = audit['cells'];mitochondrial = np.zeros(n, dtype=np.uint64)
            with np.load(reference) as expected:
                for start in range(0, n, 1024):
                    stop = min(start + 1024, n)
                    block = data.X[offset + start:offset + stop].astype(np.uint64)
                    block.sum_duplicates();block.eliminate_zeros()
                    assert np.array_equal(np.asarray(block.sum(axis=1)).ravel(), expected['totalCounts'][start:stop])
                    assert np.array_equal(np.diff(block.indptr), expected['detectedFeatures'][start:stop])
                    mitochondrial[start:stop] = np.asarray(block[:, mito].sum(axis=1)).ravel()
                path = output / (accession + '.npz')
                np.savez_compressed(path, totalCounts=expected['totalCounts'],
                                    detectedFeatures=expected['detectedFeatures'], mitochondrialCounts=mitochondrial)
            author = data.obs.n_mito_umis.iloc[offset:offset + n].to_numpy(dtype=np.uint64)
            records.append({'accession': accession, 'start': offset, 'end': offset + n,
                            'referenceSHA256': digest(path), 'sourceAuditSHA256': digest(reference),
                            'authorMitoDisagreements': int(np.count_nonzero(author != mitochondrial))})
            offset += n
            print(accession, 'verified cells', offset, flush=True)
        assert offset == prepared['cells']
        result = {'status': 'passed', 'sourceSHA256': prepared['sourceSHA256'],
                  'preparerSHA256': digest(Path(__file__)), 'cells': offset,
                  'mitochondrialFeatureCount': len(mito), 'libraries': records,
                  'authorMitoDisagreements': sum(r['authorMitoDisagreements'] for r in records),
                  'nativeExecutionQualified': False}
        (output / 'receipt.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
        print(json.dumps({k: v for k, v in result.items() if k != 'libraries'}), flush=True)
    finally:
        data.file.close()


if __name__ == '__main__':
    main()
