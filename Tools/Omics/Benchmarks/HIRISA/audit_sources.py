#!/usr/bin/env python3
"""Audit every deposited sparse count and label without fitting response models.

--available-only reports an explicitly incomplete acquisition; cached per-file
checks can be reused only when both source and auditor hashes match.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sqlite3

import h5py
import numpy as np
from scipy import sparse
from acquire import digest, one


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--available-only', action='store_true')
    args = parser.parse_args()
    root = args.root
    design = json.loads((root / 'design.json').read_text())
    audit = root / 'source-audit'
    audit.mkdir(exist_ok=True)
    auditor_sha = digest(Path(__file__))
    schema = None
    axes = None
    results = []
    missing = []
    # Disk-backed global identity check; source UUIDs remain original identities.
    database = audit / 'identities.sqlite'
    connection = sqlite3.connect(database)
    connection.execute('DROP TABLE IF EXISTS identities')
    connection.execute('CREATE TABLE identities (uuid TEXT PRIMARY KEY, accession TEXT NOT NULL)')
    for sample in design['samples']:
        accession = sample['accession']
        name = one(sample, 'supplementary_file_1').rsplit('/', 1)[1]
        path = root / 'raw' / name
        if not path.exists():
            missing.append(accession)
            continue
        source_sha = digest(path)
        record_path = audit / (accession + '.json')
        with h5py.File(path, 'r') as h:
            matrix = h['matrix']
            genes, cells = map(int, matrix['shape'][:])
            assert matrix['data'].ndim == matrix['indices'].ndim == matrix['indptr'].ndim == 1
            entries = len(matrix['data'])
            assert matrix['indices'].shape == (entries,) and matrix['indptr'].shape == (cells + 1,)
            assert matrix['data'].dtype.kind == 'u' and matrix['data'].dtype.itemsize <= 4
            feature_ids = matrix['features/id'].asstr()[:]
            names = matrix['features/name'].asstr()[:]
            assert len(feature_ids) == len(set(feature_ids)) == len(names) == genes
            current_axes = list(zip(feature_ids, names))
            if axes is None:
                axes = current_axes
            assert current_axes == axes, ('feature axis mismatch', accession)
            feature_types = matrix['features/feature_type'].asstr()[:]
            assert len(feature_types) == genes and set(feature_types) == {'Gene Expression'}
            obs = matrix['observations']
            current_schema = {k: str(v.dtype) for k, v in obs.items()}
            if schema is None:
                schema = current_schema
            assert set(current_schema) == set(schema), ('observation columns differ', accession)
            assert all(v.shape == (cells,) for v in obs.values())
            uuid = obs['cell_uuid'].asstr()[:]
            assert np.array_equal(matrix['barcodes'].asstr()[:], uuid)
            assert all(uuid) and len(set(uuid)) == cells
            try:
                connection.executemany('INSERT INTO identities VALUES (?, ?)',
                                       ((str(u), accession) for u in uuid))
            except sqlite3.IntegrityError as error:
                raise AssertionError('Repeated original cell_uuid across deposited samples') from error
            connection.commit()
            assert set(obs['pool_id'].asstr()[:]) == {one(sample, 'batch pool')}
            assert set(obs['batch_id'].asstr()[:]) == {one(sample, 'batch pool')}
            if record_path.exists():
                record = json.loads(record_path.read_text())
                assert record['sourceSHA256'] == source_sha and record['auditorSHA256'] == auditor_sha
                assert digest(audit / record['referencePath']) == record['referenceSHA256']
                results.append(record)
                print(json.dumps({'accession': accession, 'cached': True, 'cells': cells}), flush=True)
                continue
            offsets = matrix['indptr'][:].astype(np.int64)
            assert offsets[0] == 0 and offsets[-1] == entries and np.all(offsets[1:] >= offsets[:-1])
            cell_counts = np.zeros(cells, dtype=np.uint64)
            cell_features = np.zeros(cells, dtype=np.uint32)
            gene_counts = np.zeros(genes, dtype=np.uint64)
            canonical_entries = 0
            explicit_zeros = 0
            for start in range(0, cells, 1024):
                stop = min(start + 1024, cells)
                left, right = int(offsets[start]), int(offsets[stop])
                indices = matrix['indices'][left:right].astype(np.int32)
                data = matrix['data'][left:right].astype(np.uint64)
                assert np.all(indices >= 0) and np.all(indices < genes)
                # UInt32 source maxima times the complete source entry count
                # must fit UInt64 before any sparse sum is allowed.
                assert entries * int(np.iinfo(matrix['data'].dtype).max) < 2**64
                explicit_zeros += int(np.count_nonzero(data == 0))
                block = sparse.csr_matrix((data, indices, offsets[start:stop + 1] - left),
                                          shape=(stop - start, genes))
                block.sum_duplicates()
                block.eliminate_zeros()
                canonical_entries += block.nnz
                cell_counts[start:stop] = np.asarray(block.sum(axis=1)).ravel()
                cell_features[start:stop] = np.diff(block.indptr)
                gene_counts += np.asarray(block.sum(axis=0)).ravel()
            assert int(cell_counts.sum()) == int(gene_counts.sum())
            reference = audit / (accession + '.npz')
            np.savez_compressed(reference, totalCounts=cell_counts, detectedFeatures=cell_features,
                                geneCounts=gene_counts)
            labels = {}
            for key in ['celltype.l1', 'celltype.l2', 'celltype.l3']:
                labels[key] = dict(sorted(Counter(obs[key].asstr()[:]).items()))
            record = {'accession': accession, 'sourceSHA256': source_sha, 'auditorSHA256': auditor_sha,
                      'sourceBytes': path.stat().st_size, 'cells': cells, 'features': genes,
                      'storedEntries': entries, 'canonicalNonzeros': canonical_entries,
                      'explicitZeros': explicit_zeros, 'totalUMIs': int(cell_counts.sum()),
                      'countDtype': str(matrix['data'].dtype), 'observationSchema': current_schema,
                      'authorQCDisagreements': {
                          'n_umis': int(np.count_nonzero(obs['n_umis'][:] != cell_counts)),
                          'n_genes': int(np.count_nonzero(obs['n_genes'][:] != cell_features))},
                      'authorLabels': labels, 'referencePath': reference.name,
                      'referenceSHA256': digest(reference)}
            record_path.write_text(json.dumps(record, sort_keys=True, indent=2) + '\n')
            results.append(record)
            print(json.dumps({k: record[k] for k in ['accession', 'cells', 'storedEntries', 'totalUMIs',
                                                   'authorQCDisagreements']}), flush=True)
    unique_cells = connection.execute('SELECT COUNT(*) FROM identities').fetchone()[0]
    connection.close()
    summary = {'complete': len(results) == 131, 'auditedFiles': len(results), 'missingAccessions': missing,
               'auditorSHA256': auditor_sha, 'designSHA256': digest(root / 'design.json'),
               'uniqueOriginalCellUUIDs': unique_cells,
               **{k: sum(r[k] for r in results) for k in
                  ['cells', 'sourceBytes', 'storedEntries', 'canonicalNonzeros', 'totalUMIs']},
               'features': len(axes or []), 'fitsPerformed': False,
               'nativeExecutionQualified': False,
               'authorQCDisagreements': {k: sum(r['authorQCDisagreements'][k] for r in results)
                                         for k in ['n_umis', 'n_genes']}}
    assert unique_cells == summary['cells']
    output = audit / ('summary.json' if summary['complete'] else 'partial-summary.json')
    output.write_text(json.dumps(summary, sort_keys=True, indent=2) + '\n')
    print(json.dumps(summary), flush=True)
    assert args.available_only or summary['complete'], 'Complete source acquisition required'

if __name__ == '__main__':
    main()
