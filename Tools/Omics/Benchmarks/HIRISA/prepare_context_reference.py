#!/usr/bin/env python3
"""Verify frozen roles and independently aggregate sparse preparation-transfer counts."""
import argparse
from collections import Counter
import json
from pathlib import Path
import time
import resource

import h5py
import numpy as np
from scipy import sparse

from acquire import digest
from freeze_context_transfer import LINEAGES, SOURCE, encoded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    frozen = root/'context-transfer'
    receipt = json.loads((frozen/'freeze.json').read_text())
    for name, identity in receipt['files'].items():
        assert digest(frozen/name) == identity['sha256']
    assert digest(root/'hirisa.h5ad') == receipt['sourceSHA256'] == SOURCE
    membership = np.fromfile(frozen/'membership.bin', dtype=np.int8)
    aggregates = json.loads((frozen/'aggregates.json').read_text())
    folds = json.loads((frozen/'folds.json').read_text())
    keys = sorted(aggregates)
    key_index = {k: i for i, k in enumerate(keys)}
    lookup = {}
    for key, aggregate in aggregates.items():
        for accession in aggregate['accessions']:
            pair = (accession, aggregate['authorLabel'])
            assert pair not in lookup
            lookup[pair] = key_index[key]
    started = time.monotonic()
    with h5py.File(root/'hirisa.h5ad', 'r') as h:
        accessions = h['obs/geo_accession'].asstr()[:]
        labels = h['obs/celltype.l1'].asstr()[:]
        group = np.asarray([lookup.get(pair, -1) for pair in zip(accessions, labels)], dtype=np.int32)
        expected_membership = np.full(len(group), -1, dtype=np.int8)
        group_cells = np.bincount(group[group >= 0], minlength=len(keys))
        for number, (lineage, _, _) in enumerate(LINEAGES):
            plan = json.loads((frozen/(lineage+'-plan.json')).read_text())
            rows = np.flatnonzero(np.isin(group, [key_index[k] for k in keys if aggregates[k]['lineage'] == lineage]))
            assert plan['cellSelection']['observationIndices'] == rows.tolist()
            expected_membership[rows] = number
            sample_map = {s['id']: s for s in plan['mapping']['samples']}
            for key in keys:
                a = aggregates[key]
                if a['lineage'] != lineage:
                    continue
                assert group_cells[key_index[key]] == a['cells'] >= 10
                for accession in a['accessions']:
                    sample = sample_map[accession]
                    assert sample['donorID'] == a['donor']
                    assert sample['biologicalReplicateID'] == a['donor']+'-'+a['preparation']
                    assert sample['condition'] == ('none' if a['role'] == 'control' else 'IFNa')
        np.testing.assert_array_equal(membership, expected_membership)
        for fold in folds:
            training = [aggregates[k] for k in fold['trainingAggregates']]
            query, target = (aggregates[fold[k]] for k in ('queryControl', 'scoringTarget'))
            assert query['donor'] == target['donor'] == fold['heldOutDonor']
            assert query['role'] == 'control' and target['role'] == 'treated'
            assert query['preparation'] == target['preparation'] == fold['queryPreparation']
            assert query['batchPools'] == target['batchPools']
            assert all(a['donor'] != query['donor'] and a['preparation'] == fold['trainingPreparation'] for a in training)
            assert Counter(a['donor'] for a in training) == {d: 2 for d in {a['donor'] for a in training}}
            assert len(training) == 8
            for donor in {a['donor'] for a in training}:
                pair = [a for a in training if a['donor'] == donor]
                assert {a['role'] for a in pair} == {'control', 'treated'}
                assert pair[0]['batchPools'] == pair[1]['batchPools']
            train_ids = set().union(*(set(a['accessions']) for a in training))
            assert train_ids.isdisjoint(query['accessions']) and train_ids.isdisjoint(target['accessions'])
        matrix = h['X']
        assert matrix.attrs['encoding-type'] == 'csr_matrix'
        cells, features = map(int, matrix.attrs['shape'])
        assert (cells, features) == (len(group), 18082)
        ptr = matrix['indptr'][:]
        counts = np.zeros((len(keys), features), dtype=np.int64)
        totals = np.zeros(len(keys), dtype=np.int64)
        entries = blocks = 0
        for start in range(0, cells, 2048):
            stop = min(start+2048, cells)
            local = group[start:stop]
            keep = np.flatnonzero(local >= 0)
            if not len(keep):
                continue
            first, last = int(ptr[start]), int(ptr[stop])
            values = matrix['data'][first:last].astype(np.int64)
            indices = matrix['indices'][first:last]
            offset = ptr[start:stop+1]-first
            block = sparse.csr_matrix((values, indices, offset), shape=(stop-start, features))
            assignment = sparse.csr_matrix((np.ones(len(keep), dtype=np.int64),
                                            (local[keep], keep)), shape=(len(keys), stop-start))
            counts += (assignment @ block).toarray()
            # A separate integer row-total accumulation verifies denominators.
            row_total = np.asarray(block.sum(axis=1)).ravel()
            np.add.at(totals, local[keep], row_total[keep])
            entries += int(np.diff(offset)[keep].sum())
            blocks += 1
        np.testing.assert_array_equal(counts.sum(axis=1), totals)
        assert (totals > 0).all() and counts.min() >= 0
        feature_group = h['var/id']
        assert not feature_group['mask'][:].any()
        ids = feature_group['values'].asstr()[:].tolist()
        assert len(ids) == features and len(set(ids)) == features
    args.out.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(args.out/'counts.npz', counts=counts, totals=totals, cells=group_cells)
    (args.out/'identities.json').write_bytes(encoded(dict(groups=keys, features=ids)))
    status = dict(status='passed', freezeSHA256=digest(frozen/'freeze.json'),
                  sourceSHA256=SOURCE, sourceCells=cells, selectedCells=int((group >= 0).sum()),
                  groups=len(keys), features=features, selectedUMIs=int(totals.sum()),
                  selectedStoredEntries=entries, sparseBlocks=blocks,
                  all120FoldDonorAndAccessionExclusionsVerified=True,
                  originalPlanMembershipAndMetadataExact=True,
                  countsSHA256=digest(args.out/'counts.npz'), identitiesSHA256=digest(args.out/'identities.json'),
                  scriptSHA256=digest(Path(__file__)), seconds=time.monotonic()-started,
                  maximumRSSBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  responseFittingPerformed=False, heldOutResponseScoringPerformed=False)
    (args.out/'check.json').write_bytes(encoded(status))
    print(json.dumps(status))


if __name__ == '__main__':
    main()
