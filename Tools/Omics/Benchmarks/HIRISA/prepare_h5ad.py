#!/usr/bin/env python3
"""Repackage all audited 10x CSC libraries as one chunked AnnData CSR source.

CSC genes by cells has exactly the CSR storage of its cells by genes transpose.
No counts are rounded, filtered, normalized, or densified. Source files remain
unchanged and all original per-cell annotations are carried into obs.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

from acquire import digest, one


def array(group, name, shape, dtype, *, string=False):
    dataset = group.create_dataset(name, shape=shape,
                                   dtype=h5py.string_dtype('utf-8') if string else dtype,
                                   chunks=True, compression='gzip', compression_opts=1,
                                   shuffle=not string)
    dataset.attrs['encoding-type'] = 'string-array' if string else 'array'
    dataset.attrs['encoding-version'] = '0.2.0'
    return dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    design = json.loads((root / 'design.json').read_text())
    receipt = json.loads((root / 'source-receipt.json').read_text())
    audit = json.loads((root / 'source-audit/summary.json').read_text())
    assert receipt['complete'] and audit['complete'] and audit['auditedFiles'] == 131
    assert receipt['designSHA256'] == audit['designSHA256'] == digest(root / 'design.json')
    assert receipt['protocolSHA256'] == digest(Path(__file__).with_name('PROTOCOL.md'))
    rows = {s['accession']: s for s in design['samples']}
    cells, genes, entries = audit['cells'], audit['features'], audit['storedEntries']
    dtypes = {}
    count_dtype = np.dtype('uint16')
    feature_columns = None
    library_dimensions = []
    for source in receipt['files']:
        path = root / source['path']
        assert digest(path) == source['sha256']
        with h5py.File(path, 'r') as h:
            matrix = h['matrix']
            n_genes, n_cells = map(int, matrix['shape'][:])
            count_dtype = np.result_type(count_dtype, matrix['data'].dtype)
            assert n_genes == genes
            library_dimensions.append((n_cells, len(matrix['data'])))
            for key, value in matrix['observations'].items():
                dtypes[key] = np.result_type(dtypes.get(key, value.dtype), value.dtype)
            columns = {k: d.asstr()[:] for k, d in matrix['features'].items()
                       if isinstance(d, h5py.Dataset) and d.shape == (genes,)}
            if feature_columns is None:
                feature_columns = columns
            assert set(columns) == set(feature_columns)
            assert all(np.array_equal(v, feature_columns[k]) for k, v in columns.items())
    assert sum(d[0] for d in library_dimensions) == cells
    assert sum(d[1] for d in library_dimensions) == entries
    assert all(d.kind in {'S', 'u', 'i', 'f'} for d in dtypes.values())
    extras = {'geo_accession': 'accession', 'geo_donor': 'subject id',
              'geo_enrichment': 'cell type', 'geo_treatment': 'treatment',
              'geo_experiment_batch': 'batch id', 'geo_pool': 'pool id',
              'geo_batch_pool': 'batch pool'}
    assert not set(extras) & set(dtypes)
    for key, source_key in extras.items():
        values = [s['accession'] if source_key == 'accession' else one(s, source_key)
                  for s in design['samples']]
        dtypes[key] = np.dtype('S' + str(max(len(v.encode()) for v in values)))
    output = root / 'hirisa.h5ad'
    assert not output.exists(), 'Refuse overwrite of a prepared source'
    with h5py.File(output, 'x') as h:
        h.attrs['encoding-type'] = 'anndata'
        h.attrs['encoding-version'] = '0.1.0'
        x = h.create_group('X')
        x.attrs['encoding-type'] = 'csr_matrix'
        x.attrs['encoding-version'] = '0.1.0'
        x.attrs['shape'] = np.array([cells, genes], dtype=np.int64)
        data = x.create_dataset('data', shape=(entries,), dtype=count_dtype, chunks=(65536,),
                                compression='gzip', compression_opts=1, shuffle=True)
        indices = x.create_dataset('indices', shape=(entries,), dtype='int32', chunks=(65536,),
                                   compression='gzip', compression_opts=1, shuffle=True)
        pointers = x.create_dataset('indptr', shape=(cells + 1,), dtype='int64', chunks=(65536,),
                                    compression='gzip', compression_opts=1, shuffle=True)
        pointers[0] = 0
        obs = h.create_group('obs')
        obs.attrs['encoding-type'] = 'dataframe'
        obs.attrs['encoding-version'] = '0.2.0'
        obs.attrs['_index'] = '_index'
        obs.attrs['column-order'] = np.array(sorted(dtypes), dtype=h5py.string_dtype())
        array(obs, '_index', (cells,), dtypes['cell_uuid'], string=True)
        for key, dtype in dtypes.items():
            array(obs, key, (cells,), dtype, string=dtype.kind == 'S')
        frame = pd.DataFrame(feature_columns, index=pd.Index(feature_columns['id'], name='_index'))
        write_elem(h, 'var', frame)
        for group in ['obsm', 'obsp', 'varm', 'varp', 'layers']:
            write_elem(h, group, {})
        write_elem(h, 'uns', {
            'source': 'GSE306664 complete deposited labeled 10x Flex libraries',
            'source_receipt_sha256': digest(root / 'source-receipt.json'),
            'metadata_design_sha256': digest(root / 'design.json'),
            'protocol_sha256': receipt['protocolSHA256'],
            'original_feature_group_preserved_in': 'unchanged per-accession H5 source files',
            'count_transform': 'transpose CSC genes-by-cells storage into CSR cells-by-genes; integer widening only',
            'author_cell_labels_authoritative': False})
        cell_cursor = entry_cursor = 0
        for source, (n_cells, n_entries) in zip(receipt['files'], library_dimensions):
            sample = rows[source['accession']]
            with h5py.File(root / source['path'], 'r') as f:
                matrix = f['matrix']
                sl = slice(cell_cursor, cell_cursor + n_cells)
                for key in matrix['observations']:
                    values = matrix['observations'][key][:]
                    obs[key][sl] = values
                    assert np.array_equal(obs[key][sl], values, equal_nan=values.dtype.kind == 'f')
                obs['_index'][sl] = matrix['observations/cell_uuid'][:]
                for key, source_key in extras.items():
                    value = sample['accession'] if source_key == 'accession' else one(sample, source_key)
                    obs[key][sl] = value.encode()
                source_pointers = matrix['indptr'][:].astype(np.int64)
                pointers[cell_cursor + 1:cell_cursor + n_cells + 1] = source_pointers[1:] + entry_cursor
                for start in range(0, n_entries, 65536):
                    stop = min(start + 65536, n_entries)
                    destination = slice(entry_cursor + start, entry_cursor + stop)
                    values = matrix['data'][start:stop]
                    feature_indices = matrix['indices'][start:stop]
                    data[destination] = values
                    indices[destination] = feature_indices
                    assert np.array_equal(data[destination], values)
                    assert np.array_equal(indices[destination], feature_indices)
                cell_cursor += n_cells
                entry_cursor += n_entries
            print(json.dumps({'accession': source['accession'], 'cellsWritten': cell_cursor,
                              'entriesWritten': entry_cursor}), flush=True)
        assert cell_cursor == cells and entry_cursor == entries
        assert int(pointers[-1]) == entries
        # AnnData decodes both frame metadata and source IDs; matrix remains on disk.
        assert read_elem(h['var']).index.equals(frame.index)
    samples = []
    for sample in design['samples']:
        donor = one(sample, 'subject id')
        samples.append({'id': sample['accession'], 'donorID': donor,
                        'biologicalReplicateID': donor,
                        'batchID': one(sample, 'batch pool'),
                        'condition': one(sample, 'treatment'), 'organism': 'NCBITaxon:9606'})
    mapping = {'schemaVersion': 1, 'id': 'hirisa-complete-deposited-release',
               'evidence': 'measured', 'countUnit': 'umiCount', 'matrixPath': 'X',
               'sourceDescription': 'Complete 131-library GSE306664 deposited 10x Flex counts; author labels retained as predictions. Group by original accession to keep library aggregates separate; donor identities remain original.',
               'samples': samples, 'sampleColumn': 'geo_accession',
               'groupColumn': 'geo_accession', 'featureNameColumn': 'name',
               'mitochondrialFeatureIDs': [g for g, n in zip(feature_columns['id'], feature_columns['name']) if n.startswith('MT-')]}
    (root / 'stream-plan.json').write_text(json.dumps({'schemaVersion': 1, 'mapping': mapping,
                                                      'contrasts': []}, sort_keys=True, indent=2) + '\n')
    prepared = {'complete': True, 'sourceSHA256': digest(output), 'sourceBytes': output.stat().st_size,
                'sourceReceiptSHA256': digest(root / 'source-receipt.json'),
                'preparerSHA256': digest(Path(__file__)), 'cells': cells, 'features': genes,
                'storedEntries': entries, 'countValuesAndIndicesChecked': entries,
                'nativeExecutionQualified': False, 'sourceFilesPreserved': 131}
    (root / 'prepared-receipt.json').write_text(json.dumps(prepared, sort_keys=True, indent=2) + '\n')
    print(json.dumps(prepared), flush=True)

if __name__ == '__main__':
    main()
