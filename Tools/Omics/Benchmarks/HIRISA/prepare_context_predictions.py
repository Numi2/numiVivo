#!/usr/bin/env python3
"""Bind frozen context-transfer folds to verified native source groups."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import ijson


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False).encode()


def fingerprint(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def write(path, value):
    with path.open('xb') as stream:
        stream.write(canonical(value)+b'\n')


def projection(bulk, indices, lineage):
    groups, columns, counts, offsets = [], [], [], [0]
    matrix = bulk['matrix']
    for index in indices:
        groups.append({**bulk['groups'][index], 'cellGroup': lineage})
        first, last = matrix['rowOffsets'][index:index+2]
        columns.extend(matrix['featureIndices'][first:last])
        counts.extend(matrix['counts'][first:last])
        offsets.append(len(counts))
    return dict(method='explicit-source-group-projection-v1', countUnit=bulk['countUnit'],
                groups=groups, featureIDs=bulk['featureIDs'],
                matrix=dict(cellCount=len(indices), featureCount=len(bulk['featureIDs']),
                            rowOffsets=offsets, featureIndices=columns, counts=counts))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--lineage', required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    root, lineage = args.root, args.lineage
    assert not args.out.exists()
    frozen = read(root/'freeze.json')
    for name in ('folds.json', 'aggregates.json', lineage+'-plan.json'):
        assert sha(root/name) == frozen['files'][name]['sha256']
    native = root/(lineage+'-native')
    complete = read(root/(lineage+'-complete.json'))
    check = read(root/(lineage+'-independent.json'))
    assert complete['status'] == check['status'] == 'passed' and complete['nativeReplayPassed']
    assert sha(root/(lineage+'-independent.json')) == complete['independentCheckSHA256']
    assert sha(native/'report.json') == complete['reportSHA256'] == check['reportSHA256']
    assert sha(native/'receipt.json') == complete['receiptSHA256']
    assert read(native/'plan.json') == read(root/(lineage+'-plan.json'))
    aggregates = read(root/'aggregates.json')
    with (native/'report.json').open('rb') as stream:
        bulk = next(ijson.items(stream, 'pseudobulk'))
    by_key = {}
    for i, group in enumerate(bulk['groups']):
        donor, preparation = group['biologicalReplicateID'].rsplit('-', 1)
        role = {'none': 'control', 'IFNa': 'treated'}[group['condition']]
        key = f'{lineage}:{preparation}:{donor}:{role}'
        assert key not in by_key
        by_key[key] = i
        expected = aggregates[key]
        assert group['donorID'] == donor and group['cellGroup'] == expected['authorLabel']
        assert group['sampleIDs'] == expected['accessions'] and group['batchIDs'] == expected['batchPools']
    assert set(by_key) == {k for k, a in aggregates.items() if a['lineage'] == lineage}
    original_folds = [f for f in read(root/'folds.json') if f['lineage'] == lineage]
    plans, checks, scoring = [], [], []
    for fold in original_folds:
        assert fold['metadataEligible'] and not fold['unavailableTrainingAggregates']
        training = [by_key[k] for k in fold['trainingAggregates']]
        query = [by_key[fold['queryControl']]]
        target = by_key[fold['scoringTarget']]
        assert target not in training+query and len(training) == 8
        assert all(bulk['groups'][i]['donorID'] != fold['heldOutDonor'] for i in training)
        source_train = set(j for i in training for j in bulk['groups'][i]['sourceCellIndices'])
        assert source_train.isdisjoint(bulk['groups'][query[0]]['sourceCellIndices'])
        assert source_train.isdisjoint(bulk['groups'][target]['sourceCellIndices'])
        expected = [fingerprint(projection(bulk, rows, lineage)) for rows in (training, query)]
        # Perturb every stored count in every excluded aggregate, including the
        # held-out treated row, in a temporary in-memory copy only. The native
        # source and held-out outcomes are unchanged. Later compare these exact
        # projected input fingerprints to the native model/prediction receipts.
        mutant = {**bulk, 'matrix': {**bulk['matrix'], 'counts': bulk['matrix']['counts'].copy()}}
        changed = 0
        for i in range(len(bulk['groups'])):
            if i in training+query:
                continue
            first, last = bulk['matrix']['rowOffsets'][i:i+2]
            mutant['matrix']['counts'][first:last] = [x+1009 for x in bulk['matrix']['counts'][first:last]]
            changed += last-first
        assert changed > 0
        assert expected == [fingerprint(projection(mutant, rows, lineage)) for rows in (training, query)]
        plans.append(dict(id=fold['id'], perturbationID='hirisa-IFNa', controlCondition='none',
                          treatmentCondition='IFNa', cellGroup=lineage,
                          trainingGroupIndices=training, queryGroupIndices=query))
        checks.append(dict(id=fold['id'], trainingAggregateSHA256=expected[0], queryAggregateSHA256=expected[1],
                           excludedStoredCountsMutated=changed, mutatedProjectionsExact=True))
        scoring.append(dict(id=fold['id'], nativeTargetGroup=target, queryControl=fold['queryControl'], scoringTarget=fold['scoringTarget']))
    assert len(plans) == 20
    args.out.mkdir(parents=True)
    write(args.out/'plan.json', dict(schemaVersion=1, sourceReport={'bytes': list(bytes.fromhex(complete['reportSHA256']))},
          featureNamespace='ensembl-gene-id', provenance='Frozen HIRISA preparation-transfer protocol; exact author lineage, donor-excluded training and control-only queries. No held-out treated inputs. '+sha(root/'freeze.json'), folds=plans))
    write(args.out/'projection-checks.json', checks)
    write(args.out/'scoring-targets.json', scoring)
    write(args.out/'folds.json', original_folds)
    write(args.out/'transport.json', dict(status='passed', lineage=lineage, folds=20,
          frozenMetadataSHA256=sha(root/'freeze.json'), sourceReportSHA256=complete['reportSHA256'],
          nativeCompleteSHA256=sha(root/(lineage+'-complete.json')), scriptSHA256=sha(Path(__file__)),
          sourceGroupKeys=by_key, predictionFittingStarted=False,
          files={p.name: sha(p) for p in sorted(args.out.iterdir())}))
    print(json.dumps(dict(lineage=lineage, folds=20, status='passed', excludedCountMutations=20)))


if __name__ == '__main__':
    main()
