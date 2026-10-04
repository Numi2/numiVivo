"""Numerical intervention decisions, shared by the campaign and workspace.

Utilities are molecular research objectives, never efficacy or safety scores.
Cell dispersion is deliberately not used as confidence in a population mean.
"""
import numpy as np
from wetlab import require


def objective(features, specification):
    require(set(specification) == {'genes', 'preserveGenes', 'penalty'}, 'Typed molecular objective required')
    require(isinstance(specification['genes'], list) and specification['genes'], 'Select an RNA program')
    require(isinstance(specification['preserveGenes'], list), 'Preservation program required')
    require(len(set(specification['genes'])) == len(specification['genes']), 'Repeated objective gene')
    require(len(set(specification['preserveGenes'])) == len(specification['preserveGenes']), 'Repeated preservation gene')
    require(set(specification['genes'] + specification['preserveGenes']) <= set(features), 'Unmodeled objective gene')
    penalty = specification['penalty']
    require(type(penalty) in (int, float) and np.isfinite(penalty) and 0 <= penalty <= 10, 'Penalty outside [0,10]')
    return ([features.index(g) for g in specification['genes']],
            [features.index(g) for g in specification['preserveGenes']], float(penalty))


def utility(direct_delta, receiver_delta, indices):
    """Higher is better: reduce selected RNA, penalize absolute receiving effect."""
    genes, preserve, penalty = indices
    direct = -np.asarray(direct_delta)[..., genes].mean(axis=-1)
    if preserve:
        require(receiver_delta is not None, 'Receiving-population prediction unavailable; cannot evaluate preservation')
        outside = np.abs(np.asarray(receiver_delta)[..., preserve]).mean(axis=-1)
    else:
        outside = np.zeros_like(direct)
    return direct - penalty * outside, direct, outside


def selection_score(candidates, predicted, observed, training_mean):
    """Exact uniform-random expectation; candidates include the no-change option.

    The regret reference is the best *measured eligible* candidate, not an oracle
    for interventions absent from this experiment. Ties are lexical, declared.
    """
    require(len(candidates) >= 2 and len(set(candidates)) == len(candidates), 'Candidate set')
    require('no-intervention' in candidates, 'Explicit no-intervention option required')
    p, y, b = [np.asarray(v, float) for v in (predicted, observed, training_mean)]
    require(p.shape == y.shape == b.shape == (len(candidates),), 'Candidate score axes differ')
    require(all(np.isfinite(v).all() for v in (p, y, b)), 'Unavailable candidate score')
    choose = lambda v: sorted(range(len(candidates)), key=lambda i: (-v[i], candidates[i]))[0]
    selected, baseline = choose(p), choose(b)
    best, none = choose(y), candidates.index('no-intervention')
    return {'selected': candidates[selected], 'selectedObservedUtility': float(y[selected]),
            'bestObservedEligible': candidates[best], 'regret': float(y[best] - y[selected]),
            'randomExpectedUtility': float(y.mean()), 'gainOverRandom': float(y[selected] - y.mean()),
            'trainingMeanSelected': candidates[baseline], 'gainOverTrainingMean': float(y[selected] - y[baseline]),
            'gainOverNoIntervention': float(y[selected] - y[none]),
            'uncertainty': 'UNQUALIFIED: descriptive cohort result; no independent-unit interval',
            'reliableWinner': False}


def calibrated_controls(effects):
    """Development-only instrument diagnostic, never retrospective promotion."""
    y = np.asarray(effects, float)
    require(y.ndim == 2 and len(y) >= 3 and y.shape[1] >= 2, 'Calibration needs three interventions')
    controls = {'identity-positive': y, 'wrong-identity': np.roll(y, 1, axis=0),
                'wrong-magnitude': 2 * y, 'generic-average': np.tile(y.mean(0), (len(y), 1)),
                'uninformative-zero': np.zeros_like(y), 'opposite-response': -y}
    result = {}
    for name, pred in controls.items():
        correct = np.abs(y) > .1
        result[name] = {'responseRMSE': float(np.sqrt(np.mean((pred-y)**2))),
                        'directionAgreement': float(np.mean(np.sign(pred[correct]) == np.sign(y[correct]))) if correct.any() else None,
                        'meanAbsoluteMagnitudeError': float(np.mean(np.abs(np.abs(pred)-np.abs(y)))),
                        'meanProgramRegret': float(np.mean(np.max(-y, axis=0) + y[np.argmax(-pred, axis=0), np.arange(y.shape[1])]))}
    return {'status': 'passed' if all(v['responseRMSE'] > 0 for k,v in result.items() if k != 'identity-positive') else 'uninformative-development-data',
            'controls': result, 'scope': 'Metric sensitivity on exposed development responses; identity-positive is an oracle diagnostic, not a learned prediction',
            'changesV03Verdict': False}


def acquisition_backtest(predictions, truth, coverage, budget=3):
    """Retrospective acquisition: query then select the best observed candidate.

    Disagreement is a heuristic, not information gain. Uniform random's expected
    regret is calculated over fixed seeded simulations with the same budget.
    """
    p, y = np.asarray(predictions), np.asarray(truth)
    require(p.ndim == 2 and p.shape[1] == len(y), 'Acquisition candidate axis')
    budget = min(budget, len(y)); rng = np.random.default_rng(271828)
    orders = {'disagreement': np.argsort(-p.std(0), kind='stable')[:budget],
              'coverage': np.argsort(np.asarray(coverage), kind='stable')[:budget]}
    out = {k: {'selectedIndices': ix.tolist(), 'regret': float(y.max()-y[ix].max())} for k,ix in orders.items()}
    regrets = [float(y.max()-y[rng.choice(len(y), budget, replace=False)].max()) for _ in range(1000)]
    return {'budget': budget, 'policies': out, 'randomMeanRegret': float(np.mean(regrets)),
            'randomMonteCarloSD': float(np.std(regrets)), 'activeLearningBenefitEstablished': False,
            'scope': 'retrospective measurement prioritization; no model refitting, no calibrated information gain'}
