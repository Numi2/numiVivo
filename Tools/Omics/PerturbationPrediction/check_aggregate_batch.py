#!/usr/bin/env python3
"""Exercise the product batch boundary with independent synthetic count/model checks."""
import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False))


def fp(path):
    return {'bytes': list(hashlib.sha256(path.read_bytes()).digest())}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--source-binary', type=Path, required=True, help='Previously qualified product; exercises reconstruction across builds')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args();root = a.out;root.mkdir(parents=True, exist_ok=False);checks = []
    def run(name, args, success=True, binary=None):
        result = subprocess.run([str(binary or a.binary), *map(str, args)], capture_output=True, text=True)
        (root / (name + '.log')).write_text(result.stdout + result.stderr)
        assert result.returncode == (0 if success else 65), (name, result.returncode, result.stderr)
        checks.append(name)
    samples = [];obs = [];counts = []
    for donor in range(4):
        for treated in range(2):
            sample = f'd{donor}-' + ('treated' if treated else 'control')
            samples.append({'id': sample, 'donorID': f'd{donor}', 'biologicalReplicateID': f'd{donor}',
                            'condition': 'treated' if treated else 'control', 'batchID': 'synthetic', 'organism': 'NCBITaxon:9606'})
            for cell in range(4):
                obs.append({'sample': sample})
                counts.append([5 + (cell + 1) * (g + 2) + (donor + 1) ** 2 * (g % 3 + 1) + treated * (g + 1) * 3 for g in range(8)])
    counts = np.asarray(counts, dtype=np.uint64)
    data = ad.AnnData(X=csr_matrix(counts), obs=pd.DataFrame(obs, index=[f'cell-{i:02d}' for i in range(len(obs))]),
                      var=pd.DataFrame(index=[f'g{i}' for i in range(8)]))
    data.write_h5ad(root / 'source.h5ad', compression='gzip')
    mapping = {'schemaVersion': 1, 'id': 'synthetic-batch', 'evidence': 'synthetic', 'sourceDescription': 'Synthetic aggregate prediction regression',
               'countUnit': 'umiCount', 'matrixPath': 'X', 'samples': samples, 'sampleColumn': 'sample', 'groupColumn': 'sample'}
    write(root / 'source-plan.json', {'schemaVersion': 1, 'mapping': mapping, 'contrasts': []})
    source_args = ['singlecell-h5ad-pseudobulk', root / 'source.h5ad', '--plan', root / 'source-plan.json', '--output', root / 'source']
    run('source', source_args, binary=a.source_binary)
    report = json.loads((root / 'source/report.json').read_text());bulk = report['pseudobulk'];groups = bulk['groups']
    rows = {g['sampleIDs'][0]: i for i, g in enumerate(groups)}
    y = counts.reshape(8, 4, 8).sum(axis=1)
    source_y = np.array([y[[s['id'] for s in samples].index(g['sampleIDs'][0])] for g in groups])
    for row, g in enumerate(groups):
        indices = g['sourceCellIndices']
        assert np.array_equal(counts[indices].sum(axis=0), source_y[row])
        left, right = bulk['matrix']['rowOffsets'][row:row + 2]
        dense = np.zeros(8, dtype=np.uint64)
        dense[bulk['matrix']['featureIndices'][left:right]] = bulk['matrix']['counts'][left:right]
        assert np.array_equal(dense, source_y[row])
    folds = []
    for held_out in [0, 3]:
        folds.append({'id': f'hold-d{held_out}', 'perturbationID': 'synthetic-treatment', 'controlCondition': 'control',
                      'treatmentCondition': 'treated', 'cellGroup': 'synthetic-population',
                      'trainingGroupIndices': [rows[f'd{d}-{c}'] for d in range(4) if d != held_out for c in ['control', 'treated']],
                      'queryGroupIndices': [rows[f'd{held_out}-control']]})
    bad = {**folds[0], 'id': 'invalid-overlap', 'trainingGroupIndices': folds[0]['trainingGroupIndices'] + [rows['d0-treated']]}
    plan = {'schemaVersion': 1, 'sourceReport': fp(root / 'source/report.json'), 'featureNamespace': 'synthetic',
            'provenance': 'Synthetic boundary and numerical checks only', 'folds': folds + [bad]}
    write(root / 'batch-plan.json', plan)
    run('batch', ['singlecell-perturbation-batch', root / 'source', '--plan', root / 'batch-plan.json', '--output', root / 'batch'])
    receipt = json.loads((root / 'batch/receipt.json').read_text())
    assert [r['status'] for r in receipt['folds']] == ['completed', 'completed', 'failed']
    assert 'trainingAggregate' not in receipt['folds'][2]
    assert receipt['implementation'] != json.loads((root / 'source/receipt.json').read_text())['implementation']
    for i, fold in enumerate(folds):
        directory = root / 'batch/folds' / f'{i:03d}'
        model = json.loads((directory / 'model.json').read_text())
        prediction = json.loads((directory / 'prediction.json').read_text())['predictions'][0]
        training_rows = fold['trainingGroupIndices'];query_row = fold['queryGroupIndices'][0]
        yy = source_y[training_rows].astype(float);logs = np.log1p(yy / yy.sum(axis=1, keepdims=True) * 1e6)
        control, treated = logs[::2], logs[1::2];response = treated - control
        paired = yy[::2] + yy[1::2]
        selected = np.flatnonzero((paired.sum(axis=0) >= 10) & ((paired > 0).sum(axis=0) >= 2))
        centers = control[:, selected].mean(axis=0);scales = control[:, selected].std(axis=0)
        constant = np.all(control[:, selected] == control[0, selected], axis=0)
        centers[constant] = control[0, selected[constant]];scales[constant | (scales == 0)] = 1
        contexts = (control[:, selected] - centers) / scales / np.sqrt(len(selected))
        mean = response.mean(axis=0);median = np.median(response, axis=0)
        dual = np.linalg.solve(contexts @ contexts.T + np.eye(len(control)), response - mean)
        for key, expected in [('contextCenters', centers), ('contextScales', scales), ('contexts', contexts),
                              ('meanResponse', mean), ('medianResponse', median), ('dualCoefficients', dual)]:
            assert np.allclose(model[key], expected, rtol=1e-11, atol=1e-11), (i, key)
        assert model['selectedFeatureIndices'] == selected.tolist()
        assert model['source'] == receipt['folds'][i]['trainingAggregate']
        assert model['trainingDonors'] == sorted({groups[j]['donorID'] for j in training_rows})
        q = source_y[query_row].astype(float);q = np.log1p(q / q.sum() * 1e6)
        context = (q[selected] - centers) / scales / np.sqrt(len(selected))
        ridge = mean + (context @ contexts.T) @ dual
        expected_responses = {'noChange': np.zeros(8), 'meanResponse': mean, 'medianResponse': median, 'contextRidge': ridge}
        assert np.allclose(prediction['control'], q, rtol=1e-12, atol=1e-12)
        assert prediction['group']['sourceCellIndices'] == groups[query_row]['sourceCellIndices']
        for estimate in prediction['estimates']:
            response = expected_responses[estimate['baseline']];predicted = np.maximum(0, q + response)
            assert np.allclose(estimate['unclippedResponse'], response, rtol=1e-11, atol=1e-11)
            assert np.allclose(estimate['predictedTreated'], predicted, rtol=1e-11, atol=1e-11)
            assert np.allclose(estimate['predictedResponse'], predicted - q, rtol=1e-11, atol=1e-11)
            assert np.isclose(estimate['impliedCPMSum'], np.expm1(predicted).sum(), rtol=1e-11)
    checks.append('independent-count-membership-and-model-oracle')
    run('verify', ['singlecell-perturbation-batch-verify', root / 'batch'])
    run('repeat', ['singlecell-perturbation-batch', root / 'source', '--plan', root / 'batch-plan.json', '--output', root / 'repeat'])
    for path in (root / 'batch/folds').rglob('*.json'):
        assert path.read_bytes() == (root / 'repeat/folds' / path.relative_to(root / 'batch/folds')).read_bytes()
    checks.append('byte-exact-repeat')
    run('overwrite', ['singlecell-perturbation-batch', root / 'source', '--plan', root / 'batch-plan.json', '--output', root / 'batch'], False)
    forged = root / 'forged-model';shutil.copytree(root / 'batch', forged)
    path = forged / 'folds/000/model.json';model = json.loads(path.read_text());model['meanResponse'][0] += 1;write(path, model)
    outer = json.loads((forged / 'receipt.json').read_text());outer['folds'][0]['model'] = fp(path)
    write(forged / 'receipt.json', outer);write(forged / 'folds/000/receipt.json', outer['folds'][0])
    run('reject-rehashed-model', ['singlecell-perturbation-batch-verify', forged], False)
    forged = root / 'forged-source';shutil.copytree(root / 'source', forged)
    path = forged / 'report.json';changed = json.loads(path.read_text());changed['pseudobulk']['matrix']['counts'][0] += 1;write(path, changed)
    old = json.loads((forged / 'receipt.json').read_text());old['report'] = fp(path);write(forged / 'receipt.json', old)
    write(root / 'forged-plan.json', {**plan, 'sourceReport': fp(path)})
    run('reject-rehashed-source', ['singlecell-perturbation-batch', forged, '--plan', root / 'forged-plan.json', '--output', root / 'rejected-source'], False)
    assert not (root / 'rejected-source').exists()
    fabricated = root / 'fabricated-failed-fold';shutil.copytree(root / 'batch', fabricated)
    write(fabricated / 'folds/002/model.json', {})
    run('reject-fabricated-failed-fold', ['singlecell-perturbation-batch-verify', fabricated], False)
    result = {'status': 'passed', 'checks': checks, 'sourceReceipt': json.loads((root / 'source/receipt.json').read_text()),
              'batchReceipt': receipt, 'scope': 'Synthetic isolation, independent numerical oracle and replay/tamper gates; no biological prediction claim.'}
    write(root / 'checks.json', result)
    print(json.dumps({'status': 'passed', 'checks': len(checks)}))


if __name__ == '__main__':
    main()
