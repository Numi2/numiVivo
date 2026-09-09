#!/usr/bin/env python3
"""Qualify native condition-composition commands against all frozen Norman predictions."""
import argparse
import copy
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

import combinations as reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('binary', 'pseudobulk', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    base = Path(__file__).parent / 'evidence/2026-09-09-combinations'
    expected = reference.load(base / 'inputs/training.npz')
    queries = json.loads((base / 'inputs/queries.json').read_text())
    expected_report = 'c36a13a80371e80fc5c0eb0dee598cdaed69dc955c7df312a4174588bffcee9c'
    assert reference.sha(args.pseudobulk / 'report.json') == expected_report
    args.out.mkdir(parents=True, exist_ok=False)
    context = dict(id='Norman2019-K562-pooled', organism='NCBITaxon:9606', featureNamespace='Ensembl-gene',
                   perturbationNamespace='Norman2019-guide-target', countUnit='umiCount')
    selection = dict(schemaVersion=1, context=context, controlCondition='control',
                     targets=[dict(id=str(t), condition=str(t)) for t in expected['targetIDs']],
                     provenance='Complete Norman 2019 filtered release, Zenodo 13350497. Verified native aggregates; '
                     'one pooled K562 context, no independent donor/replication claim. Control and 105 declared singles only.')
    plan = dict(schemaVersion=1, context=context, queries=queries)
    reference.write(args.out / 'selection.json', selection)
    reference.write(args.out / 'queries.json', plan)
    commands = []

    def run(name, command, success=True):
        with (args.out / (name + '.log')).open('w') as log:
            result = subprocess.run(['/usr/bin/time', '-l', str(args.binary.resolve()), *map(str, command)], stdout=log, stderr=subprocess.STDOUT)
        commands.append(dict(name=name, returnCode=result.returncode, expectedSuccess=success))
        reference.write(args.out / 'commands.json', commands)
        assert result.returncode == (0 if success else 65), (name, result.returncode)

    training_path, model_path, prediction_path = (args.out / name for name in ('training.json', 'model', 'prediction'))
    run('prepare', ['singlecell-composition-prepare', args.pseudobulk, '--plan', args.out / 'selection.json', '--output', training_path])
    training = json.loads(training_path.read_text())
    assert training['featureIDs'] == expected['featureIDs'].tolist()
    assert training['selection']['context'] == context and training['evidence'] == 'measured'
    assert training['selection']['targets'] == selection['targets']
    assert bytes(training['sourceReport']['bytes']).hex() == expected_report
    assert bytes(training['source']['bytes']).hex() == 'efde6f5301fe256725dce1d980f37bd96a13481a9a16135515897368e631affc'
    matrix = training['matrix']
    assert matrix['cellCount'] == 106 and matrix['featureCount'] == 33694
    raw = np.vstack([expected['controlCounts'], expected['singleCounts']])
    for row in range(len(raw)):
        indices = np.flatnonzero(raw[row]);start, end = matrix['rowOffsets'][row:row + 2]
        assert matrix['featureIndices'][start:end] == indices.tolist()
        assert matrix['counts'][start:end] == raw[row, indices].tolist()
    run('fit', ['singlecell-composition-fit', training_path, '--output', model_path])
    run('model-verify', ['singlecell-composition-verify', model_path])
    run('predict', ['singlecell-composition-predict', model_path, '--plan', args.out / 'queries.json', '--output', prediction_path])
    run('prediction-verify', ['singlecell-composition-prediction-verify', prediction_path])
    native_model = json.loads((model_path / 'model.json').read_text())
    assert native_model['targetIDs'] == expected['targetIDs'].tolist()
    expected_cpm = raw.astype(float) / raw.sum(axis=1, keepdims=True) * 1e6
    native_cpm = np.array([native_model['controlCPM'], *native_model['targetCPM']])
    np.testing.assert_allclose(native_cpm, expected_cpm, rtol=1e-12, atol=1e-9)
    report = json.loads((prediction_path / 'report.json').read_text())
    assert report['featureIDs'] == expected['featureIDs'].tolist() and report['queryIDs'] == [q['id'] for q in queries]
    expected_diagnostics = json.loads((base / 'predictions/diagnostics.json').read_text())
    comparisons = []
    assert [m['method'] for m in report['methods']] == list(reference.METHODS)
    for method in report['methods']:
        name = method['method'];saved = reference.load(base / 'predictions' / (name + '.npz'))
        assert method['queryRows'] == saved['queryRows'].tolist()
        values = np.array(method['expression'])
        np.testing.assert_allclose(values, saved['expression'], rtol=1e-12, atol=1e-12)
        max_error = float(np.max(abs(values - saved['expression'])))
        cpm_error = 0.0
        clip_differences = []
        for i, diagnostic in enumerate(method['diagnostics']):
            prior = expected_diagnostics[i]['models'][name]
            assert diagnostic['renormalized'] is False
            error = abs(diagnostic['impliedCPMSum'] - prior['impliedCPMSum'])
            cpm_error = max(cpm_error, error)
            np.testing.assert_allclose(diagnostic['impliedCPMSum'], prior['impliedCPMSum'], rtol=1e-12, atol=1e-8)
            if diagnostic['clippedFeatures'] != prior['clippedFeatures']:
                clip_differences.append(dict(query=queries[i]['id'], native=diagnostic['clippedFeatures'], external=prior['clippedFeatures']))
        comparisons.append(dict(method=name, maximumPredictionError=max_error, maximumCPMSumError=cpm_error, clipDifferences=clip_differences))
    reference.write(args.out / 'numerical-comparison.json', comparisons)
    # Diagnostic classification differences require inspection, not silent tolerance.
    assert not any(c['clipDifferences'] for c in comparisons), 'clipping diagnostics differ; inspect numerical-comparison.json'
    del report, native_cpm, native_model
    reverse = copy.deepcopy(plan)
    for query in reverse['queries']:
        query['targets'].reverse()
    reference.write(args.out / 'reverse.json', reverse)
    run('reverse-targets', ['singlecell-composition-predict', model_path, '--plan', args.out / 'reverse.json', '--output', args.out / 'reverse'])
    assert reference.sha(args.out / 'reverse/report.json') == reference.sha(prediction_path / 'report.json')
    reordered = copy.deepcopy(training)
    reordered['selection']['targets'].reverse()
    rows = [0, *range(105, 0, -1)];offsets, columns, values = [0], [], []
    for row in rows:
        start, end = matrix['rowOffsets'][row:row + 2]
        columns.extend(matrix['featureIndices'][start:end]);values.extend(matrix['counts'][start:end]);offsets.append(len(values))
    reordered['matrix'].update(rowOffsets=offsets, featureIndices=columns, counts=values)
    reference.write(args.out / 'reordered-training.json', reordered)
    run('reordered-training', ['singlecell-composition-fit', args.out / 'reordered-training.json', '--output', args.out / 'reordered-model'])
    assert reference.sha(args.out / 'reordered-model/model.json') == reference.sha(model_path / 'model.json')
    del training, reordered, matrix, columns, values
    for name, field, value in [('context', 'id', 'other-context'), ('organism', 'organism', 'NCBITaxon:10090'),
                               ('feature-namespace', 'featureNamespace', 'other-features'),
                               ('target-namespace', 'perturbationNamespace', 'other-targets'), ('count-unit', 'countUnit', 'readCount')]:
        wrong = copy.deepcopy(plan);wrong['context'][field] = value
        path = args.out / (name + '.json');reference.write(path, wrong)
        run('reject-' + name, ['singlecell-composition-predict', model_path, '--plan', path, '--output', args.out / ('rejected-' + name)], False)
    for name in ('unknown-target', 'repeated-target', 'duplicate-query', 'empty-query', 'unknown-field'):
        wrong = copy.deepcopy(plan)
        if name == 'unknown-target':wrong['queries'][0]['targets'][0] = 'unobserved-target'
        elif name == 'repeated-target':wrong['queries'][0]['targets'][1] = wrong['queries'][0]['targets'][0]
        elif name == 'duplicate-query':wrong['queries'][1]['id'] = wrong['queries'][0]['id']
        elif name == 'empty-query':wrong['queries'] = []
        else:wrong['pairedOutcomes'] = [1]
        path = args.out / (name + '.json');reference.write(path, wrong)
        run('reject-' + name, ['singlecell-composition-predict', model_path, '--plan', path, '--output', args.out / ('rejected-' + name)], False)
    run('reject-model-overwrite', ['singlecell-composition-fit', training_path, '--output', model_path], False)
    run('reject-prediction-overwrite', ['singlecell-composition-predict', model_path, '--plan', args.out / 'queries.json', '--output', prediction_path], False)
    tampered = args.out / 'tampered-model';shutil.copytree(model_path, tampered)
    data = (tampered / 'model.json').read_bytes()
    data, count = re.subn(rb'("controlCPM":\[)([^,\]]+)', lambda m: m[1] + str(float(m[2]) + 1).encode(), data, count=1)
    assert count == 1
    (tampered / 'model.json').write_bytes(data)
    receipt = json.loads((tampered / 'receipt.json').read_text());receipt['model'] = dict(bytes=list(hashlib.sha256(data).digest()))
    reference.write(tampered / 'receipt.json', receipt)
    run('reject-rehashed-model', ['singlecell-composition-verify', tampered], False)
    assert not list(args.out.glob('.numivivo-composition-*'))
    result = dict(status='passed', queries=131, features=33694, methods=6, exactTrainingCounts=True,
                  trainingOrderExact=True, targetOrderExact=True, comparisons=comparisons,
                  binarySHA256=reference.sha(args.binary), sourceReportSHA256=expected_report,
                  nativeModelSHA256=reference.sha(model_path / 'model.json'), nativeReportSHA256=reference.sha(prediction_path / 'report.json'),
                  commands=commands, scope='Native single-context composition baselines. No unseen-target, genetic-interaction, donor, Bayesian or mechanistic qualification.')
    reference.write(args.out / 'checks.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
