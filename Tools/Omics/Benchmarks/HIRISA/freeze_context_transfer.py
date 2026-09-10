#!/usr/bin/env python3
"""Freeze complete HIRISA preparation-transfer membership without RNA values."""
import argparse
from collections import Counter, defaultdict
import copy
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from acquire import digest, one, parse_samples

HERE = Path(__file__).resolve().parent
SOURCE = '0873e698ebf8770a54e6dba09724ffbeda5e1a67dbf24d4e223552d4aca7969c'
DESIGN = 'ace6e6eb5dad1ff8fcc43d748753b05a9e6d9a7d989a14efc6fe54cb066885e8'
LINEAGES = [('B', 'B', 'Bcell'), ('Mono', 'Mono', 'Monocyte'), ('NK', 'NK', 'NK'),
            ('CD4-T', 'CD4 T', 'Tcell'), ('CD8-T', 'CD8 T', 'Tcell'),
            ('other-T', 'other T', 'Tcell')]


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def library_pairs(samples, design):
    enriched = {c['population']: c['pairs'] for c in design['comparisons'] if c['treatment'] == 'IFNa'}
    by_key = defaultdict(lambda: {'control': [], 'treated': []})
    for sample in samples:
        if one(sample, 'cell type') != 'PBMC':
            continue
        treatment = one(sample, 'treatment')
        if treatment not in ('culture_IFNa', 'culture_no_stim'):
            continue
        key = tuple(one(sample, k) for k in ('subject id', 'batch id', 'pool id'))
        by_key[key]['treated' if treatment == 'culture_IFNa' else 'control'].append(sample['accession'])
    pbmc = []
    for (donor, batch, pool), pair in sorted(by_key.items()):
        assert len(pair['control']) == len(pair['treated']) == 1, (donor, batch, pool, pair)
        pbmc.append(dict(donor=donor, batch=batch, pool=pool,
                         control=pair['control'][0], treated=pair['treated'][0]))
    return enriched, pbmc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    assert not args.out.exists()
    assert digest(root/'design.json') == DESIGN
    design = json.loads((root/'design.json').read_text())
    samples = parse_samples(root/'series.soft.gz')
    assert samples == design['samples']
    assert digest(root/'series.soft.gz') == design['sourceMetadataSHA256']
    prepared = json.loads((root/'prepared-receipt.json').read_text())
    assert prepared['sourceSHA256'] == SOURCE and prepared['cells'] == 1612594
    assert digest(root/'hirisa.h5ad') == SOURCE
    source = json.loads((root/'source-receipt.json').read_text())
    assert digest(root/'source-receipt.json') == prepared['sourceReceiptSHA256']
    enriched, pbmc = library_pairs(samples, design)
    donors = sorted({one(s, 'subject id') for s in samples})
    assert len(donors) == 5 and len(pbmc) == 9
    sample_by_id = {s['accession']: s for s in samples}
    with h5py.File(root/'hirisa.h5ad', 'r') as h:
        accessions = h['obs/geo_accession'].asstr()[:]
        labels = h['obs/celltype.l1'].asstr()[:]
        uuids = h['obs/cell_uuid'].asstr()[:]
        obs_donors = h['obs/geo_donor'].asstr()[:]
        obs_treatments = h['obs/geo_treatment'].asstr()[:]
        assert len(accessions) == len(labels) == len(uuids) == 1612594
        # Original library boundaries and annotations are checked directly;
        # counts and derived response scores are never opened.
        start = 0
        for entry in source['files']:
            path = root/entry['path']
            assert digest(path) == entry['sha256']
            with h5py.File(path, 'r') as raw:
                original = raw['matrix/observations']
                lab = original['celltype.l1'].asstr()[:]
                identity = original['cell_uuid'].asstr()[:]
            stop = start + len(lab)
            assert (accessions[start:stop] == entry['accession']).all()
            np.testing.assert_array_equal(labels[start:stop], lab)
            np.testing.assert_array_equal(uuids[start:stop], identity)
            sample = sample_by_id[entry['accession']]
            assert (obs_donors[start:stop] == one(sample, 'subject id')).all()
            assert (obs_treatments[start:stop] == one(sample, 'treatment')).all()
            start = stop
        assert start == len(labels)
    selected = np.full(len(labels), -1, dtype=np.int8)
    plans, aggregates, folds = {}, {}, []
    membership = []
    mapping_base = json.loads((root/'stream-plan.json').read_text())['mapping']
    for lineage_index, (lineage, label, population) in enumerate(LINEAGES):
        pairs = dict(enriched=enriched[population], PBMC=pbmc)
        roles = {}
        for preparation, records in pairs.items():
            for pair in records:
                for role in ('control', 'treated'):
                    accession = pair[role]
                    assert accession not in roles
                    roles[accession] = (pair['donor'], preparation, role)
        group_rows = defaultdict(list)
        indices = np.flatnonzero(np.isin(accessions, list(roles)) & (labels == label))
        assert (selected[indices] == -1).all()
        selected[indices] = lineage_index
        for accession, (donor, preparation, role) in sorted(roles.items()):
            rows = np.flatnonzero((accessions == accession) & (labels == label))
            key = f'{lineage}:{preparation}:{donor}:{role}'
            group_rows[key].extend(rows.tolist())
        for key, rows in sorted(group_rows.items()):
            rows.sort()
            contributing = sorted(set(accessions[rows]))
            lineage_id, preparation, donor, role = key.split(':')
            aggregates[key] = dict(id=key, lineage=lineage_id, authorLabel=label,
                                   preparation=preparation, donor=donor, role=role,
                                   accessions=contributing, cells=len(rows),
                                   batchPools=sorted({one(sample_by_id[a], 'batch pool') for a in contributing}),
                                   originalRowsSHA256=hashlib.sha256(np.asarray(rows, dtype='<u4').tobytes()).hexdigest())
        mapping = copy.deepcopy(mapping_base)
        mapping.update(id='hirisa-context-'+lineage, groupColumn='celltype.l1',
                       sourceDescription='Complete source retained; exact author-label preparation-transfer selection. See frozen CONTEXT_TRANSFER_PROTOCOL.md.')
        for sample in mapping['samples']:
            if sample['id'] in roles:
                donor, preparation, role = roles[sample['id']]
                sample['biologicalReplicateID'] = donor+'-'+preparation
                sample['condition'] = 'none' if role == 'control' else 'IFNa'
        plan = dict(schemaVersion=1, mapping=mapping, contrasts=[],
                    cellSelection=dict(source={'bytes': list(bytes.fromhex(SOURCE))},
                                       observationIndices=indices.tolist(),
                                       provenance='Frozen exact accession and author celltype.l1 membership; all matched PBMC pools pooled only within donor/preparation/condition; no expression filtering.'))
        plans[lineage] = encoded(plan)
        assert len(plans[lineage]) <= 2097152, (lineage, len(plans[lineage]), 'Native plan capacity requires an owner repair; do not reduce the cohort')
        membership.append(dict(lineage=lineage, authorLabel=label, enrichment=population,
                               cells=len(indices), matchedLibraryPairs=pairs,
                               planBytes=len(plans[lineage])))
        for destination in ('PBMC', 'enriched'):
            other = 'enriched' if destination == 'PBMC' else 'PBMC'
            for donor in donors:
                for kind, training_preparation in [('cross', other), ('within', destination)]:
                    training = []
                    unavailable = []
                    for train_donor in donors:
                        if train_donor == donor:
                            continue
                        keys = [f'{lineage}:{training_preparation}:{train_donor}:{role}' for role in ('control', 'treated')]
                        if all(aggregates[k]['cells'] >= 10 for k in keys):
                            training.extend(keys)
                        else:
                            unavailable.extend(keys)
                    query = f'{lineage}:{destination}:{donor}:control'
                    target = f'{lineage}:{destination}:{donor}:treated'
                    assert all(aggregates[k]['donor'] != donor for k in training)
                    folds.append(dict(id=f'{kind}-{other}-to-{destination}-{lineage}-{donor}' if kind == 'cross' else f'{kind}-{destination}-{lineage}-{donor}',
                                      kind=kind, lineage=lineage, heldOutDonor=donor,
                                      trainingPreparation=training_preparation, queryPreparation=destination,
                                      trainingAggregates=training, unavailableTrainingAggregates=unavailable,
                                      queryControl=query, scoringTarget=target,
                                      metadataEligible=(len(training) >= 4 and aggregates[query]['cells'] >= 10 and aggregates[target]['cells'] >= 10)))
    ledger = []
    for (accession, label), cells in sorted(Counter(zip(accessions, labels)).items()):
        mask = (accessions == accession) & (labels == label)
        assignments = Counter(selected[mask].tolist())
        assert len(assignments) == 1
        assignment = next(iter(assignments))
        ledger.append(dict(accession=accession, authorLabel=label, cells=cells,
                           selectedLineage=LINEAGES[assignment][0] if assignment >= 0 else None,
                           exclusionReason=None if assignment >= 0 else 'outside-fixed-library-and-author-lineage-mapping'))
    assert len(folds) == 120 and len({f['id'] for f in folds}) == 120
    assert sum(x['cells'] for x in ledger) == len(labels)
    assert digest(root/'hirisa.h5ad') == SOURCE
    args.out.mkdir(parents=True)
    for lineage, raw in plans.items():
        (args.out/(lineage+'-plan.json')).write_bytes(raw)
    (args.out/'membership.bin').write_bytes(selected.tobytes())
    for name, value in [('aggregates.json', aggregates), ('folds.json', folds), ('ledger.json', ledger), ('lineages.json', membership)]:
        (args.out/name).write_bytes(encoded(value))
    repo = HERE.parents[3]
    owners = ['VivoPerturbation.swift', 'VivoPerturbationIO.swift', 'VivoPerturbationAggregateBatch.swift', 'VivoH5ADPseudobulk.swift']
    result = dict(schemaVersion=1, status='frozen-metadata-native-count-and-prediction-execution-pending',
                  protocolSHA256=digest(HERE/'CONTEXT_TRANSFER_PROTOCOL.md'), freezerSHA256=digest(Path(__file__)),
                  sourceSHA256=SOURCE, designSHA256=DESIGN, sourceReceiptSHA256=digest(root/'source-receipt.json'),
                  sourceCells=len(labels), selectedCells=int((selected >= 0).sum()), excludedCells=int((selected < 0).sum()),
                  sourceLibrariesVerified=131, originalLabelsAndUUIDsExact=True,
                  crossPreparationFolds=60, withinPreparationReferences=60,
                  metadataEligibleFolds=sum(f['metadataEligible'] for f in folds),
                  expressionValuesRead=False, predictionFittingStarted=False,
                  modelOwnerSHA256={name: digest(repo/'Sources/NumiVivoKit/Omics'/name) for name in owners},
                  files={p.name: dict(bytes=p.stat().st_size, sha256=digest(p)) for p in sorted(args.out.iterdir())})
    (args.out/'freeze.json').write_bytes(encoded(result))
    print(json.dumps({k: v for k, v in result.items() if k not in ('files', 'modelOwnerSHA256')}))


if __name__ == '__main__':
    main()
