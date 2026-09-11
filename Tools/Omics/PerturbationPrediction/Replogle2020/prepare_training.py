#!/usr/bin/env python3
"""Prepare all five native training contexts; freeze count folds without fitting.

Guide IDs remain experimental labels, not verified gene descriptors. This tool
does not retrieve GO, fit a model, predict outcomes, or select favorable folds.
"""
import argparse
import copy
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np


BINARY = '36df30cd11afffb9bca0f88369af2dca3260759ddd9945c2c02d36073b1574db'
RNA = '7170656afc5fd23fb7ee749467e5cb60832cc1766100aba8c99ec1ea96bfa17d'
REPORT = '4c4142c29304d3f480bf5165232ad83d84d9acea0a8f536d604a9e4d62adec54'
CONTROLS = {'sgNegCtrl2', 'sgNegCtrl3'}
BASE = '72af0e8b02f57685605e077e180eb29df0866176'


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def write(path, value):
    path.write_bytes(encoded(value))


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1_048_576), b''):
            digest.update(block)
    return digest.hexdigest()


def subset(training, keep):
    """Same count-only exclusion contract as the qualified Norman driver."""
    matrix = training['matrix']
    offsets, indices, counts = [0], [], []
    for row in [0, *[i+1 for i in keep]]:
        start, end = matrix['rowOffsets'][row:row+2]
        indices.extend(matrix['featureIndices'][start:end])
        counts.extend(matrix['counts'][start:end])
        offsets.append(len(counts))
    return dict(training, selection=dict(training['selection'],
                targets=[training['selection']['targets'][i] for i in keep]),
                matrix=dict(matrix, cellCount=len(keep)+1, rowOffsets=offsets,
                            featureIndices=indices, counts=counts))


def dense_row(matrix, row):
    start, end = matrix['rowOffsets'][row:row+2]
    indices = matrix['featureIndices'][start:end]
    assert len(indices) == len(set(indices))
    result = np.zeros(matrix['featureCount'], dtype=np.uint64)
    result[np.asarray(indices, dtype=np.int64)] = matrix['counts'][start:end]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.out.resolve()
    binary = root/'runtime/numivivo-omics'
    assert sha(binary) == BINARY
    assert sha(root/'rna-projection/projected.h5ad') == RNA
    assert sha(root/'selected-aggregate/report.json') == REPORT
    assert json.loads((root/'selected-checks.json').read_text())['status'] == 'passed'
    out.mkdir(parents=True, exist_ok=False)
    before = os.statvfs(out)
    original = json.loads((root/'selected-aggregate/report.json').read_text())
    plan = json.loads((root/'selected-plan.json').read_text())
    assert plan['cellSelection']['observationIndices'] == original['sourceObservationIndices']
    groups = original['pseudobulk']['groups']
    gemgroups = sorted({g['cellGroup'] for g in groups})
    assert gemgroups == ['gemgroup-'+str(i) for i in range(1, 6)]
    guides = sorted({g['condition'] for g in groups} - CONTROLS)
    assert len(guides) == 30 and len(groups) == 160
    assert all({g['condition'] for g in groups if g['cellGroup'] == gem} == set(guides)|CONTROLS for gem in gemgroups)
    plan['mapping']['id'] = 'replogle2020-upr-controls-pooled-within-gemgroup'
    plan['mapping']['sourceDescription'] += '; both paper-defined controls pooled only within original gemgroup; guide identities remain unresolved at gene level'
    for sample in plan['mapping']['samples']:
        condition = 'pooled-control' if sample['condition'] in CONTROLS else sample['condition']
        sample['condition'] = sample['batchID']+'|'+condition
    write(out/'pooled-plan.json', plan)
    commands = []

    def run(label, command, rejection=None):
        start = time.monotonic()
        with (out/(label+'.log')).open('wb') as log:
            result = subprocess.run(['/usr/bin/time', '-l', str(binary), *map(str, command)],
                stdout=log, stderr=subprocess.STDOUT, env=dict(os.environ,
                NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib'))
        commands.append(dict(label=label, command=list(map(str, command)),
            returnCode=result.returncode, seconds=time.monotonic()-start, expectedRejection=rejection))
        write(out/'commands.json', commands)
        if rejection:
            assert result.returncode != 0 and rejection in (out/(label+'.log')).read_text()
        else:
            assert result.returncode == 0, label
        print(label, result.returncode, flush=True)

    bundle = out/'pooled-aggregate'
    run('aggregate', ['singlecell-h5ad-pseudobulk', root/'rna-projection/projected.h5ad',
        '--plan', out/'pooled-plan.json', '--output', bundle])
    run('aggregate-verify', ['singlecell-h5ad-pseudobulk-verify', bundle])
    report = json.loads((bundle/'report.json').read_text())
    bulk = report['pseudobulk']
    assert report['sourceObservationIndices'] == original['sourceObservationIndices']
    assert report['quality'] == original['quality']
    assert report['canonicalNonzeros'] == original['canonicalNonzeros'] == 105512809
    assert bulk['featureIDs'] == original['pseudobulk']['featureIDs']
    assert len(bulk['groups']) == 155
    reference = np.load(root/'selected-reference.npz')
    inventory = json.loads((root/'prepared-complete/inventory.json').read_text())
    reference_rows = {name: i for i, name in enumerate(inventory['groups'])}
    audit, lookup = [], {}
    for index, group in enumerate(bulk['groups']):
        gem, condition = group['condition'].split('|')
        assert group['cellGroup'] == gem and group['batchIDs'] == [gem]
        assert group['biologicalReplicateID'] == 'K562-pooled-replication-unresolved' and 'donorID' not in group
        labels = CONTROLS if condition == 'pooled-control' else {condition}
        inputs = [g for g in groups if g['cellGroup'] == gem and g['condition'] in labels]
        assert len(inputs) == len(labels)
        expected = sum((reference['counts'][reference_rows[gem+'|'+label]].astype(np.uint64)
                       for label in sorted(labels)), np.zeros(33694, dtype=np.uint64))
        np.testing.assert_array_equal(dense_row(bulk['matrix'], index), expected)
        assert group['sourceCellIndices'] == sorted(i for g in inputs for i in g['sourceCellIndices'])
        assert group['sampleIDs'] == sorted(gem+'|'+label for label in labels)
        lookup[group['condition']] = index
        audit.append(dict(gemgroup=gem, condition=condition, sourceGuideLabels=sorted(labels),
            cells=len(group['sourceCellIndices']), totalUMI=int(expected.sum()), allGeneCountsExact=True))
    write(out/'pooling-audit.json', audit)
    del original, reference
    folds, contexts = [], []
    for gem in gemgroups:
        selection = dict(schemaVersion=1, context=dict(id='Replogle2020-K562-'+gem,
            organism='NCBITaxon:9606', featureNamespace='Ensembl-gene',
            perturbationNamespace='Replogle2020-original-guide-label', countUnit='umiCount'),
            controlCondition=gem+'|pooled-control',
            targets=[dict(id=guide, condition=gem+'|'+guide) for guide in guides],
            provenance='Complete confident Replogle cohort; paper-defined controls pooled within original technical gemgroup. Guide labels are not verified gene identities. No model fit or score.')
        write(out/(gem+'-selection.json'), selection)
        training_path = out/(gem+'-training.json')
        run(gem+'-prepare', ['singlecell-composition-prepare', bundle, '--plan',
            out/(gem+'-selection.json'), '--output', training_path])
        training = json.loads(training_path.read_text())
        assert training['selection'] == selection and training['evidence'] == 'measured'
        assert bytes(training['source']['bytes']).hex() == RNA
        assert bytes(training['sourceReport']['bytes']).hex() == sha(bundle/'report.json')
        assert training['featureIDs'] == bulk['featureIDs'] and training['matrix']['cellCount'] == 31
        for row, condition in enumerate([selection['controlCondition'], *[x['condition'] for x in selection['targets']]]):
            np.testing.assert_array_equal(dense_row(training['matrix'], row), dense_row(bulk['matrix'], lookup[condition]))
        contexts.append(dict(gemgroup=gem, trainingSHA256=sha(training_path), rows=31, features=33694,
            controlCells=next(x['cells'] for x in audit if x['gemgroup']==gem and x['condition']=='pooled-control')))
        for held, guide in enumerate(guides):
            keep = [i for i in range(len(guides)) if i != held]
            selected = subset(training, keep)
            mutation = dict(training, matrix=dict(training['matrix'], counts=training['matrix']['counts'].copy()))
            lo, hi = mutation['matrix']['rowOffsets'][held+1:held+3]
            mutation['matrix']['counts'][lo:hi] = [1]*(hi-lo)
            assert subset(mutation, keep) == selected
            assert selected['matrix']['cellCount'] == 30
            assert guide not in [x['id'] for x in selected['selection']['targets']]
            folds.append(dict(gemgroup=gem, heldGuide=guide, trainingGuides=[guides[i] for i in keep],
                trainingRows=[0, *[i+1 for i in keep]], selectedTrainingSHA256=hashlib.sha256(encoded(selected)).hexdigest(),
                heldCountMutationInvariant=True, descriptorIdentityStatus='unverified', predictionStatus='not-run'))
        if gem == gemgroups[0]:
            negative = copy.deepcopy(selection)
            negative['targets'][0]['condition'] = gemgroups[1]+'|'+guides[0]
            write(out/'cross-context-rejection.json', negative)
            run('cross-context-rejection', ['singlecell-composition-prepare', bundle, '--plan',
                out/'cross-context-rejection.json', '--output', out/'must-not-exist.json'],
                'composition conditions cross biological contexts')
            assert not (out/'must-not-exist.json').exists()
    assert len(folds) == 150
    write(out/'count-folds.json', dict(status='count-selection-frozen-descriptors-pending',
        frozenUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(), contexts=contexts, folds=folds,
        predictionsFrozen=False, fitted=False, scored=False))
    after = os.statvfs(out)
    artifacts = [dict(path=str(p.relative_to(out)), bytes=p.stat().st_size, SHA256=sha(p))
        for p in sorted(out.rglob('*')) if p.is_file() and p.suffix != '.h5ad']
    write(out/'checks.json', dict(status='passed-count-preparation-only', sourceBase=BASE,
        binarySHA256=BINARY, RNA_SHA256=RNA, originalSelectedReportSHA256=REPORT,
        scriptSHA256=sha(Path(__file__)), selectedCells=32829, features=33694, originalGuideGroups=160,
        pooledGroups=155, contexts=contexts, frozenCountFolds=150, commands=len(commands),
        exactIndependentCounts=True, exactMembership=True, exactCellQC=True, nativeReplay=True,
        crossContextRejected=True, allHeldCountMutationsInvariant=True,
        descriptorIdentitiesVerified=False, predictionsFrozen=False, fitted=False, scored=False,
        freeBytesBefore=before.f_bavail*before.f_frsize, freeBytesAfter=after.f_bavail*after.f_frsize,
        files=artifacts))
    print('All five contexts and 150 count exclusions qualified; descriptors and prediction remain pending.', flush=True)


if __name__ == '__main__':
    main()
