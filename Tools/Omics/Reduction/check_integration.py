#!/usr/bin/env python3
"""Pinned Harmony reference, donor mixing, and measured response preservation.

Operates on the complete retained cohort and its existing sparse RNA/PCA report.
Condition labels are evaluation-only for Harmony; the erasure control deliberately
uses them and is never a candidate integration method.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import harmonypy
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def neighbors(x, groups=None, k=30):
    """Exact neighbors with stable ties, excluding self; no dense N x N matrix."""
    result = np.empty((len(x), k), dtype=np.int64)
    index = np.arange(len(x))
    for start in range(0, len(x), 64):
        distances = cdist(x[start:start + 64], x)
        for offset, row in enumerate(distances):
            i = start + offset
            valid = index != i
            if groups is not None:
                valid &= groups == groups[i]
            candidates = index[valid]
            if len(candidates) < k:
                raise ValueError('Insufficient cells for condition-stratified neighbors')
            result[i] = candidates[np.lexsort((candidates, row[candidates]))[:k]]
    return result


def metrics(x, donors, conditions, program):
    within = neighbors(x, conditions)
    global_neighbors = neighbors(x)
    levels = np.unique(donors)
    entropy = []
    same_donor = []
    expected_same = []
    for i, indices in enumerate(within):
        p = np.array([np.mean(donors[indices] == level) for level in levels])
        entropy.append(float(-np.sum(p[p > 0] * np.log(p[p > 0])) / np.log(len(levels))))
        same_donor.append(float(np.mean(donors[indices] == donors[i])))
        eligible = conditions == conditions[i]
        expected_same.append(float((np.sum(eligible & (donors == donors[i])) - 1) / (eligible.sum() - 1)))
    folds = []
    for donor in levels:
        test = donors == donor
        model = make_pipeline(StandardScaler(), LogisticRegression(C=1, max_iter=2000, solver='lbfgs'))
        model.fit(x[~test], conditions[~test])
        if model[-1].n_iter_.max() >= 2000:
            raise ValueError('Condition classifier did not converge')
        prediction = model.predict(x[test])
        folds.append(dict(donor=str(donor), cells=int(test.sum()),
                          balancedAccuracy=float(balanced_accuracy_score(conditions[test], prediction)),
                          solverIterations=int(model[-1].n_iter_.max())))
    # All summaries below weight donors equally; the cell counts are also retained.
    def donor_mean(values):
        return float(np.mean([np.mean(np.asarray(values)[donors == donor]) for donor in levels]))
    smoothed = program[global_neighbors].mean(axis=1)
    within_smoothed = program[within].mean(axis=1)
    program_correlations = []
    within_condition_correlations = []
    for donor in levels:
        mask = donors == donor
        program_correlations.append(dict(donor=str(donor),
            spearman=float(spearmanr(program[mask], smoothed[mask]).statistic)))
        for condition in np.unique(conditions):
            selected = mask & (conditions == condition)
            within_condition_correlations.append(dict(donor=str(donor), condition=str(condition),
                spearman=float(spearmanr(program[selected], within_smoothed[selected]).statistic)))
    return dict(conditionStratifiedDonorEntropy=donor_mean(entropy),
                conditionStratifiedSameDonorFraction=donor_mean(same_donor),
                compositionExpectedSameDonorFraction=donor_mean(expected_same),
                sameDonorExcess=donor_mean(np.asarray(same_donor) - expected_same),
                conditionNeighborPurity=donor_mean((conditions[global_neighbors] == conditions[:, None]).mean(axis=1)),
                leaveOneDonorOutBalancedAccuracy=float(np.mean([f['balancedAccuracy'] for f in folds])),
                donorFolds=folds, programNeighborSpearman=float(np.mean([f['spearman'] for f in program_correlations])),
                programDonorCorrelations=program_correlations,
                withinConditionProgramSpearman=float(np.mean([f['spearman'] for f in within_condition_correlations])),
                withinConditionProgramCorrelations=within_condition_correlations)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--program', nargs='+', default=['IFI6', 'IFIT1', 'ISG15', 'MX1', 'ISG20'])
    args = parser.parse_args()
    if importlib.metadata.version('harmonypy') != '2.0.0':
        raise ValueError('Reference requires harmonypy==2.0.0')
    raw = args.report.read_bytes()
    report = json.loads(raw)
    dataset = report['processed']['dataset']
    reduction = report['reduction']
    cells = reduction['cells']
    identity = lambda cell: (cell['sampleID'], cell['barcode'])
    if list(map(identity, cells)) != list(map(identity, dataset['cells'])):
        raise ValueError('PCA and retained RNA cell axes differ')
    samples = {sample['id']: sample for sample in dataset['samples']}
    donors = np.asarray([samples[cell['sampleID']].get('donorID', '') for cell in cells])
    conditions = np.asarray([samples[cell['sampleID']].get('condition', '') for cell in cells])
    if '' in donors or '' in conditions or len(np.unique(donors)) < 3 or len(np.unique(conditions)) != 2:
        raise ValueError('Benchmark requires at least three known donors and two conditions')
    composition = pd.crosstab(donors, conditions)
    if (composition == 0).any().any():
        raise ValueError('Every donor must have both conditions; confounded designs are ineligible')
    x = np.asarray(reduction['scores'], dtype=np.float64)
    if x.shape[0] != len(cells) or not np.isfinite(x).all():
        raise ValueError('Invalid PCA matrix')
    normalized = report['processed']['normalized']
    matrix = sparse.csr_matrix((normalized['values'], normalized['featureIndices'], normalized['rowOffsets']),
                              shape=(len(cells), normalized['featureCount']))
    feature_ids = [feature['id'] for feature in dataset['features']]
    if len(set(args.program)) != len(args.program) or any(feature_ids.count(gene) != 1 for gene in args.program):
        raise ValueError('Program requires unique, present feature IDs; no silent marker dropping')
    program = np.asarray(matrix[:, [feature_ids.index(gene) for gene in args.program]].mean(axis=1)).ravel()
    baseline = metrics(x, donors, conditions, program)
    # Label-informed negative control must reveal loss even if mixing looks good.
    erased = x.copy()
    for condition in np.unique(conditions):
        mask = conditions == condition
        erased[mask] -= erased[mask].mean(axis=0)
    control = metrics(erased, donors, conditions, program)
    options = dict(theta=2, lamb=1, sigma=0.1, nclust=min(round(len(x) / 30), 100),
                   tau=0, block_size=0.05, max_iter_harmony=10, max_iter_kmeans=4,
                   epsilon_cluster=0.001, epsilon_harmony=0.01, alpha=0.2,
                   batch_prop_cutoff=1e-5, ncores=1)
    runs = []
    def gates(measurement):
        return dict(mixingImproved=measurement['sameDonorExcess'] < baseline['sameDonorExcess'],
                    conditionPreserved=measurement['leaveOneDonorOutBalancedAccuracy'] >= baseline['leaveOneDonorOutBalancedAccuracy'] - 0.02,
                    programPreserved=measurement['programNeighborSpearman'] >= baseline['programNeighborSpearman'] - 0.05,
                    withinConditionProgramPreserved=measurement['withinConditionProgramSpearman'] >= baseline['withinConditionProgramSpearman'] - 0.05)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    for seed in [7, 19, 41]:
        result = harmonypy.run_harmony(x, pd.DataFrame({'donor': donors}), ['donor'],
                                      **options, random_state=seed, verbose=False)
        corrected = np.asarray(result.Z_corr)
        if corrected.shape != x.shape or not np.isfinite(corrected).all():
            raise ValueError('Invalid reference corrected coordinates')
        scores_path = args.out.with_name(args.out.stem + f'-seed{seed}.npz')
        np.savez_compressed(scores_path, scores=corrected)
        measurement = metrics(corrected, donors, conditions, program)
        # Predeclared engineering margins for this one cohort, not biology proof.
        runs.append(dict(seed=seed, metrics=measurement, gates=gates(measurement),
                         objectiveHarmony=list(result.objective_harmony), objectiveClustering=list(result.objective_kmeans),
                         clusteringRounds=list(result.kmeans_rounds), coordinateFile=scores_path.name,
                         coordinateSHA256=hashlib.sha256(scores_path.read_bytes()).hexdigest()))
    summary = dict(schemaVersion=1, status='reference-baseline-measured',
                   reportSHA256=hashlib.sha256(raw).hexdigest(), cells=len(cells), dimensions=x.shape[1],
                   composition={str(d): {str(c): int(composition.loc[d, c]) for c in composition.columns} for d in composition.index},
                   options=options, evaluationNeighbors=30,
                   gateMargins=dict(maximumConditionAccuracyLoss=0.02, maximumProgramCorrelationLoss=0.05,
                                    minimumErasureControlAccuracyLoss=0.10),
                   program=dict(featureIDs=args.program, definition='Unweighted mean of retained log-normalized RNA; no fitted weights'),
                   baseline=baseline, conditionErasureControl=control, references=runs,
                   controlDetectsLoss=control['leaveOneDonorOutBalancedAccuracy'] < baseline['leaveOneDonorOutBalancedAccuracy'] - 0.1,
                   qualification='One paired B-cell cohort. Transductive integration sees all unlabeled cells; donor-held-out classifier is not prospective perturbation prediction. No native integration qualification or multi-cell-type preservation claim.',
                   versions={name: importlib.metadata.version(name) for name in ['harmonypy', 'numpy', 'scipy', 'pandas', 'scikit-learn']})
    if 'integration' in report:
        native = report['integration']
        if native['cells'] != cells or native['options']['covariate'] != 'donor':
            raise ValueError('Native integration axes or covariate differ')
        corrected = np.asarray(native['scores'])
        if corrected.shape != x.shape or not np.isfinite(corrected).all():
            raise ValueError('Invalid native corrected coordinates')
        measurement = metrics(corrected, donors, conditions, program)
        np.testing.assert_allclose(native['relativeImprovements'],
            -np.diff(native['objectives']) / np.maximum(np.abs(native['objectives'][:-1]), 1e-12), rtol=1e-12, atol=1e-12)
        if native['stoppingReason'] == 'relative-objective-tolerance':
            assert 0 <= native['relativeImprovements'][-1] < native['options']['relativeTolerance']
        elif native['stoppingReason'] == 'objective-increase':
            assert native['relativeImprovements'][-1] < 0
        else:
            assert native['stoppingReason'] == 'iteration-limit'
            assert len(native['relativeImprovements']) == native['options']['maximumIterations']
        # Independent dense (B+1)-square solve; never a cells-by-genes matrix.
        memberships = np.asarray(native['memberships'])
        batch = np.asarray(native['cellLevels'])
        assert [native['levels'][i] for i in batch] == donors.tolist()
        np.testing.assert_allclose(memberships.sum(axis=1), 1, rtol=0, atol=1e-12)
        assert memberships.shape == (len(x), native['options']['clusters'])
        assert np.isfinite(memberships).all() and (memberships >= 0).all()
        reconstructed = x.copy()
        observed = np.stack([memberships[batch == b].sum(axis=0) for b in range(len(native['levels']))], axis=1)
        sizes = np.bincount(batch)
        expected = memberships.sum(axis=0)[:, None] * sizes[None, :] / len(x)
        distances = 2 * (1 - np.asarray(native['assignmentScores']) @ np.asarray(native['assignmentCenters']).T)
        entropy = np.sum(memberships[memberships > 0] * np.log(memberships[memberships > 0]))
        objective = (np.sum(memberships * distances) + native['options']['temperature'] * entropy +
                     native['options']['temperature'] * native['options']['diversity'] *
                     np.sum(observed * np.log((observed + expected + 1) / (2 * expected + 1)))) * 2000 / len(x)
        np.testing.assert_allclose(objective, native['objectives'][-1], rtol=1e-10, atol=1e-10)
        for cluster in range(memberships.shape[1]):
            active = np.flatnonzero(observed[cluster] / sizes > 1e-5)
            if len(active) < 2:
                continue
            masses = observed[cluster, active]
            sums = np.stack([memberships[batch == b, cluster] @ x[batch == b] for b in active])
            normal = np.diag(np.r_[masses.sum(), masses + native['options']['ridge']])
            normal[0, 1:] = masses
            normal[1:, 0] = masses
            beta = np.linalg.solve(normal, np.vstack([sums.sum(axis=0), sums]))
            for slot, b in enumerate(active):
                selected = batch == b
                reconstructed[selected] -= memberships[selected, cluster, None] * beta[slot + 1]
        np.testing.assert_allclose(reconstructed, corrected, rtol=1e-10, atol=1e-10)
        summary['native'] = dict(metrics=measurement, gates=gates(measurement),
                                stoppingReason=native['stoppingReason'], objectives=native['objectives'],
                                maximumRidgeResidual=native['maximumRidgeResidual'],
                                independentObjective=float(objective),
                                independentCorrectionMaximumError=float(np.max(np.abs(reconstructed - corrected))))
        summary['status'] = 'native-and-reference-measured'
        summary['qualification'] = 'One paired B-cell cohort; independent native initialization, not coordinate identity with Harmony. Transductive integration sees all unlabeled cells. Donor-held-out classifier is not prospective perturbation prediction or multi-cell-type preservation proof.'
    args.out.write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: summary[key] for key in ['status', 'cells', 'controlDetectsLoss']}))
    for run in runs:
        print(json.dumps(dict(seed=run['seed'], gates=run['gates'], metrics={key: value for key, value in run['metrics'].items() if not isinstance(value, list)})))
    if 'native' in summary:
        print(json.dumps({key: value for key, value in summary['native'].items() if key != 'metrics'}))
        if not all(summary['native']['gates'].values()):
            raise ValueError('Native preservation/mixing gate failed; retain evidence')
    if not summary['controlDetectsLoss']:
        raise ValueError('Biological preservation evaluation failed its loss control')


if __name__ == '__main__':
    main()
