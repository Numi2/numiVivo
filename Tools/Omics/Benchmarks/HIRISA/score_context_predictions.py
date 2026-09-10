#!/usr/bin/env python3
"""Independently reconstruct frozen native context predictions, then score all folds."""
import argparse
from pathlib import Path
import json
import numpy as np
from prepare_context_predictions import sha, read, write

METHODS = ('noChange', 'meanResponse', 'medianResponse', 'contextRidge')
LINEAGES = ('B', 'Mono', 'NK', 'CD4-T', 'CD8-T', 'other-T')


def reconstruct(model, query, training_counts, control_counts):
    y = training_counts.astype(float)
    logs = np.log1p(y/y.sum(axis=1, keepdims=True)*1e6)
    controls, response = logs[::2], logs[1::2]-logs[::2]
    paired = y[::2]+y[1::2]
    selected = np.flatnonzero((paired.sum(axis=0) >= 10) & ((paired > 0).sum(axis=0) >= 2))
    assert selected.tolist() == model['selectedFeatureIndices']
    centers = controls[:, selected].mean(axis=0)
    scales = controls[:, selected].std(axis=0)
    constant = np.all(controls[:, selected] == controls[0, selected], axis=0)
    centers[constant] = controls[0, selected[constant]]
    scales[constant | (scales == 0)] = 1
    contexts = (controls[:, selected]-centers)/scales/np.sqrt(len(selected))
    mean, median = response.mean(axis=0), np.median(response, axis=0)
    dual = np.linalg.solve(contexts@contexts.T+np.eye(len(controls)), response-mean)
    errors = {}

    def close(name, actual, expected):
        actual, expected = np.asarray(actual), np.asarray(expected)
        assert actual.shape == expected.shape and np.isfinite(actual).all(), name
        errors[name] = float(np.max(np.abs(actual-expected), initial=0))
        assert np.allclose(actual, expected, rtol=1e-9, atol=1e-9), (name, errors[name])

    for name, value in [('contextCenters', centers), ('contextScales', scales), ('contexts', contexts),
                        ('meanResponse', mean), ('medianResponse', median), ('dualCoefficients', dual)]:
        close(name, model[name], value)
    q_counts = control_counts.astype(float)
    assert query['libraryCounts'] == int(q_counts.sum())
    q = np.log1p(q_counts/q_counts.sum()*1e6)
    close('control', query['control'], q)
    context = (q[selected]-centers)/scales/np.sqrt(len(selected))
    expected = dict(noChange=np.zeros(len(q)), meanResponse=mean, medianResponse=median,
                    contextRidge=mean+(context@contexts.T)@dual)
    estimates = {e['baseline']: e for e in query['estimates']}
    assert set(estimates) == set(METHODS)
    for name, effect in expected.items():
        e = estimates[name]
        treated = np.maximum(0, q+effect)
        close(name+'-unclippedResponse', e['unclippedResponse'], effect)
        close(name+'-treated', e['predictedTreated'], treated)
        close(name+'-response', e['predictedResponse'], treated-q)
        close(name+'-impliedCPMSum', e['impliedCPMSum'], np.expm1(treated).sum())
    return q, selected, estimates, errors


def metrics(query, target, estimates, selected):
    truth = target-query
    result = []
    for family, indices in [('all-source-genes', np.arange(len(query))), ('training-selected-context-genes', selected)]:
        for method in METHODS:
            estimate = estimates[method]
            treated = np.asarray(estimate['predictedTreated'])[indices]
            response = np.asarray(estimate['predictedResponse'])[indices]
            error, delta = treated-target[indices], response-truth[indices]
            pearson = float(np.corrcoef(response, truth[indices])[0, 1]) if np.ptp(response) > 0 and np.ptp(truth[indices]) > 0 else None
            assert pearson is None or np.isfinite(pearson)
            result.append(dict(method=method, family=family, features=len(indices),
                treatedRMSE=float(np.sqrt(np.mean(error**2))), treatedMAE=float(np.mean(np.abs(error))),
                responseRMSE=float(np.sqrt(np.mean(delta**2))), responseMAE=float(np.mean(np.abs(delta))),
                responsePearson=pearson, impliedCPMSum=estimate['impliedCPMSum'],
                clippedGenes=int(np.count_nonzero(query+np.asarray(estimate['unclippedResponse']) < 0))))
    return result


def summarize(results):
    summaries, comparisons = [], []
    for lineage in LINEAGES:
        for destination in ('PBMC', 'enriched'):
            subset = [r for r in results if r['lineage'] == lineage and r['queryPreparation'] == destination]
            for kind in ('cross', 'within'):
                rows = [r for r in subset if r['kind'] == kind]
                assert len(rows) == 5
                complete = all(r['status'] == 'completed' for r in rows)
                for family in ('all-source-genes', 'training-selected-context-genes'):
                    for method in METHODS:
                        scores = [s for r in rows for s in r['scores'] if s['family'] == family and s['method'] == method]
                        means = {}
                        for name in ('treatedRMSE', 'treatedMAE', 'responseRMSE', 'responseMAE', 'responsePearson', 'impliedCPMSum'):
                            values = [s[name] for s in scores]
                            means[name] = float(np.mean(values)) if complete and all(v is not None for v in values) else None
                        summaries.append(dict(lineage=lineage, queryPreparation=destination, kind=kind,
                            family=family, method=method, expectedFolds=5, completedFolds=len(scores),
                            complete=complete, equalDonorMeans=means))
            rows = [s for s in summaries if s['lineage'] == lineage and s['queryPreparation'] == destination and s['family'] == 'all-source-genes']
            values = {(s['kind'], s['method']): s['equalDonorMeans']['responseRMSE'] for s in rows}
            complete = all(s['complete'] for s in rows)
            comparisons.append(dict(lineage=lineage, queryPreparation=destination, complete=complete,
                ridgeBeatsNoChangeAndCrossMean=(values['cross', 'contextRidge'] < values['cross', 'noChange'] and values['cross', 'contextRidge'] < values['cross', 'meanResponse']) if complete else None,
                crossMinusWithinResponseRMSE={m: values['cross', m]-values['within', m] if complete else None for m in METHODS},
                comparisonFamily='all-source-genes; secondary training-selected panels differ across training preparations'))
    return summaries, comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    freeze = read(root/'prediction-output-freeze.json')
    assert freeze['targetsOpenedForScoring'] is False and freeze['folds'] == 120
    assert freeze['frozenMetadataSHA256'] == sha(root/'freeze.json')
    for name, digest in freeze['files'].items():
        assert sha(root/name) == digest, name
    # The complete output freeze is verified before any scoring count array is opened.
    reference = read(root/'reference/check.json')
    assert reference['freezeSHA256'] == sha(root/'freeze.json')
    assert sha(root/'reference/counts.npz') == reference['countsSHA256']
    assert sha(root/'reference/identities.json') == reference['identitiesSHA256']
    identities = read(root/'reference/identities.json')
    counts = np.load(root/'reference/counts.npz', allow_pickle=False)['counts']
    group_index = {k: i for i, k in enumerate(identities['groups'])}
    all_folds = read(root/'folds.json')
    assert len(all_folds) == 120 and sha(root/'folds.json') == read(root/'freeze.json')['files']['folds.json']['sha256']
    results = []
    for lineage in LINEAGES:
        bundle = root/(lineage+'-prediction')
        inputs = root/'prediction-inputs'/lineage
        transport = read(inputs/'transport.json')
        for name, digest in transport['files'].items():
            assert sha(inputs/name) == digest
        source_folds = [f for f in all_folds if f['lineage'] == lineage]
        assert read(inputs/'folds.json') == source_folds
        plan, receipt = read(bundle/'plan.json'), read(bundle/'receipt.json')
        assert plan == read(inputs/'plan.json')
        assert receipt['plan']['bytes'] == list(bytes.fromhex(sha(bundle/'plan.json')))
        assert receipt['sourceReport'] == plan['sourceReport'] == {'bytes': list(bytes.fromhex(transport['sourceReportSHA256']))}
        assert receipt['source']['bytes'] == list(bytes.fromhex(read(root/'freeze.json')['sourceSHA256']))
        projected = read(inputs/'projection-checks.json')
        assert len(receipt['folds']) == len(projected) == len(source_folds) == 20
        for i, (fold, status, check) in enumerate(zip(source_folds, receipt['folds'], projected)):
            assert fold['id'] == status['id'] == check['id'] == plan['folds'][i]['id']
            directory = bundle/'folds'/f'{i:03d}'
            assert read(directory/'receipt.json') == status
            base = {k: fold[k] for k in ('id', 'lineage', 'kind', 'heldOutDonor', 'trainingPreparation', 'queryPreparation')}
            if status['status'] != 'completed':
                assert status['status'] == 'failed' and status.get('failure')
                results.append(dict(**base, status='failed', failure=status['failure'], scores=[]))
                continue
            model, prediction = read(directory/'model.json'), read(directory/'prediction.json')
            for name in ('model', 'prediction'):
                assert status[name]['bytes'] == list(bytes.fromhex(sha(directory/(name+'.json'))))
            assert model['method'] == 'paired-donor-log1p-CPM-response-baselines-alpha1-v1'
            assert model['featureIDs'] == prediction['featureIDs'] == identities['features']
            assert model['source'] == status['trainingAggregate'] == prediction['referenceSource']
            assert status['trainingAggregate']['bytes'] == list(bytes.fromhex(check['trainingAggregateSHA256']))
            assert status['queryAggregate']['bytes'] == list(bytes.fromhex(check['queryAggregateSHA256']))
            assert check['mutatedProjectionsExact']
            donors = sorted({k.split(':')[2] for k in fold['trainingAggregates']})
            assert model['trainingDonors'] == donors and fold['heldOutDonor'] not in donors
            expected_order = [f'{lineage}:{fold["trainingPreparation"]}:{donor}:{role}' for donor in donors for role in ('control', 'treated')]
            assert expected_order == fold['trainingAggregates']
            assert len(prediction['predictions']) == 1
            query = prediction['predictions'][0]
            assert query['group']['donorID'] == fold['heldOutDonor'] and query['group']['cellGroup'] == lineage
            q, selected, estimates, errors = reconstruct(model, query,
                counts[[group_index[k] for k in expected_order]], counts[group_index[fold['queryControl']]])
            # Select the held-out treated aggregate only after this fold's
            # implementation and excluded-input fingerprints have been checked.
            observed = counts[group_index[fold['scoringTarget']]].astype(float)
            assert observed.sum() > 0
            target = np.log1p(observed/observed.sum()*1e6)
            results.append(dict(**base, status='completed', scores=metrics(q, target, estimates, selected),
                                numericalMaximumAbsoluteErrors=errors, excludedInputFingerprintExact=True))
        print(lineage, 'checked', flush=True)
    summaries, comparisons = summarize(results)
    for name, digest in freeze['files'].items():
        assert sha(root/name) == digest, 'Frozen output changed during scoring'
    write(args.out, dict(schemaVersion=1, status='passed' if all(r['status'] == 'completed' for r in results) else 'incomplete',
        folds=results, summaries=summaries, comparisons=comparisons, foldCount=len(results),
        completedFolds=sum(r['status'] == 'completed' for r in results), outputFreezeSHA256=sha(root/'prediction-output-freeze.json'),
        referenceCheckSHA256=sha(root/'reference/check.json'), scorerSHA256=sha(Path(__file__)),
        scope='Frozen cross-preparation and matched within-preparation donor-excluded RNA response predictions; annotation and batch confounding retained; no general biological qualification'))
    print(json.dumps(dict(folds=len(results), contrasts=len(comparisons), ridgePrimaryWins=sum(c['ridgeBeatsNoChangeAndCrossMean'] is True for c in comparisons))))


if __name__ == '__main__':
    main()
