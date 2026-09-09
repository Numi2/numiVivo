#!/usr/bin/env python3
"""Audit the complete scPerturb Norman filtered release, without cell-matrix densification."""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import h5py
import numpy as np

SOURCE = {
    'url': 'https://zenodo.org/records/13350497/files/NormanWeissman2019_filtered.h5ad',
    'bytes': 698680199,
    'md5': 'c870e6967d91c017d9da827bab183cd6',
    'sha256': 'efde6f5301fe256725dce1d980f37bd96a13481a9a16135515897368e631affc',
}


def column(frame, name):
    value = frame[name]
    if isinstance(value, h5py.Group):
        assert value.attrs['encoding-type'] == 'categorical'
        labels = value['categories'].asstr()[:]
        codes = value['codes'][:]
        assert np.all((codes >= 0) & (codes < len(labels)))
        return labels[codes]
    return value.asstr()[:] if h5py.check_string_dtype(value.dtype) else value[:]


def frequencies(values):
    labels, counts = np.unique(values, return_counts=True)
    return {str(k): int(v) for k, v in zip(labels, counts)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    md5, sha = hashlib.md5(), hashlib.sha256()
    with args.source.open('rb') as stream:
        for block in iter(lambda: stream.read(8388608), b''):
            md5.update(block)
            sha.update(block)
    assert args.source.stat().st_size == SOURCE['bytes']
    assert md5.hexdigest() == SOURCE['md5'] and sha.hexdigest() == SOURCE['sha256']
    with h5py.File(args.source, 'r') as f:
        x, obs, var = f['X'], f['obs'], f['var']
        assert x.attrs['encoding-type'] == 'csc_matrix'
        assert list(x.attrs['shape']) == [111445, 33694]
        barcodes = column(obs, '_index').astype(str)
        genes = column(var, 'ensemble_id').astype(str)
        symbols = column(var, '_index').astype(str)
        for axis in (barcodes, genes, symbols):
            assert len(set(axis)) == len(axis)
        perturbations = column(obs, 'perturbation').astype(str)
        conditions, codes = np.unique(perturbations, return_inverse=True)
        guides = column(obs, 'guide_id')
        control_guides = frequencies(guides[perturbations == 'control'])
        assert all(all(target.startswith('NegCtrl') for target in guide.split(';')[0].split('_'))
                   for guide in control_guides)
        assert not any('no_reads_found' in guide for guide in guides)
        singles = sorted(c for c in conditions if c != 'control' and '_' not in c)
        pairs = sorted(c for c in conditions if '_' in c)
        assert len(singles) == 105 and len(pairs) == 131 and len(conditions) == 237
        assert all(len(c.split('_')) == 2 and set(c.split('_')) <= set(singles) for c in pairs)
        offsets = x['indptr'][:]
        assert len(offsets) == len(genes) + 1 and offsets[0] == 0
        assert np.all(offsets[1:] >= offsets[:-1])
        assert offsets[-1] == len(x['data']) == len(x['indices']) == 361582621
        # Dense condition x gene reference only: no dense cell x gene matrix.
        counts = np.zeros((len(conditions), len(genes)), dtype=np.uint64)
        totals = np.zeros(len(barcodes), dtype=np.uint64)
        detected = np.zeros(len(barcodes), dtype=np.uint64)
        mitochondrial = np.zeros(len(barcodes), dtype=np.uint64)
        mito = np.char.startswith(symbols, 'MT-')
        maximum = 0
        for j in range(len(genes)):
            start, end = int(offsets[j]), int(offsets[j + 1])
            rows, values = x['indices'][start:end], x['data'][start:end]
            assert np.all((rows >= 0) & (rows < len(barcodes)))
            # This pinned release is sorted with no duplicate coordinates/zeros.
            assert np.all(rows[1:] > rows[:-1])
            assert np.all(np.isfinite(values) & (values > 0) & (values == np.floor(values)))
            integer = values.astype(np.uint64)
            maximum = max(maximum, int(integer.max(initial=0)))
            np.add.at(counts[:, j], codes[rows], integer)
            totals[rows] += integer
            detected[rows] += 1
            if mito[j]:
                mitochondrial[rows] += integer
        assert int(counts.sum()) == int(totals.sum()) == 1635387239
        np.savez_compressed(args.out / 'reference.npz', counts=counts, conditions=conditions,
                            feature_ids=genes, gene_symbols=symbols, barcodes=barcodes,
                            cell_conditions=perturbations, totals=totals, detected=detected,
                            mitochondrial=mitochondrial)
        mapping = dict(schemaVersion=1, id='norman2019-full-filtered', evidence='measured',
                       sourceDescription='Norman 2019 doi:10.1126/science.aax4438; complete scPerturb filtered release, Zenodo 13350497. '
                       'K562 CRISPRa pooled condition aggregates. Independent biological replication and donor identity are unresolved; '
                       'eight GEM groups are not biological replicates. No differential-expression inference requested.',
                       countUnit='umiCount', matrixPath='X', sampleColumn='perturbation',
                       groupColumn='cell_line', featureIDColumn='ensemble_id',
                       mitochondrialFeatureIDs=genes[mito].tolist(),
                       samples=[dict(id=str(c), condition=str(c), biologicalReplicateID='K562-pooled-replication-unresolved',
                                     batchID='pooled-eight-gemgroups', organism='NCBITaxon:9606') for c in conditions])
        (args.out / 'plan.json').write_text(json.dumps(dict(schemaVersion=1, mapping=mapping, contrasts=[]), indent=2) + '\n')
        splits = dict(schemaVersion=1, source=SOURCE, control='control',
                      unseenCombinations=dict(trainingConditions=singles, queryConditions=pairs,
                                              rule='All pairs held out together; both target singles observed. Control is available.'),
                      unseenTargets=[dict(target=t,
                                          trainingConditions=[str(c) for c in conditions if c != 'control' and t not in c.split('_')],
                                          queryConditions=[str(c) for c in conditions if t in c.split('_')]) for t in singles],
                      scope='Condition-mean prediction; neither cell nor GEM-group splits establish donor generalization. '
                      'These are split definitions, not fitted or validated predictions.')
        (args.out / 'splits.json').write_text(json.dumps(splits, indent=2) + '\n')
        audit = dict(source=SOURCE, shape=[len(barcodes), len(genes)], storedEntries=int(offsets[-1]),
                     aggregateNonzeros=int(np.count_nonzero(counts)), totalCounts=int(totals.sum()),
                     maximumCount=maximum, uniqueBarcodes=True, uniqueFeatureIDs=True,
                     sortedCoordinates=True, duplicateCoordinates=0, nonpositiveCounts=0,
                     nonfiniteCounts=0, fractionalCounts=0, mitochondrialFeatureCount=int(mito.sum()),
                     conditions=frequencies(perturbations), controlGuides=control_guides,
                     gemgroups=frequencies(column(obs, 'gemgroup')), singles=len(singles), pairs=len(pairs),
                     python=platform.python_version(), numpy=np.__version__, h5py=h5py.__version__,
                     seconds=time.monotonic() - started)
        (args.out / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
        print(json.dumps(dict(status='passed', cells=len(barcodes), genes=len(genes), conditions=len(conditions))))


if __name__ == '__main__':
    main()
