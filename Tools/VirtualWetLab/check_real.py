#!/usr/bin/env python3
"""Bounded real-data workflow acceptance; independent NumPy oracle, all donors.

This is software and development-data evidence, not fresh biological validation.
Run once against an empty destination. No tuning or success-dependent retries.
"""
import argparse
from pathlib import Path
import numpy as np
import anndata as ad
from wetlab import catalog, predict, reveal, verify, summary, read, write, load_assay

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--assay', type=Path, required=True); p.add_argument('--binary', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=False)
cat = catalog(a.assay)
write(a.output / 'acceptance-plan.json', {'scope': 'workflow acceptance on previously inspected Kang data',
    'specimens': [s['donor'] for s in cat['specimens']], 'stoppingRule': 'All eight donors, once each; no tuning',
    'numericalTolerance': 'absolute 1e-9, relative 1e-9',
    'checks': ['exact aggregation from original measured counts', 'independent ridge solve',
               'all baseline predictions', 'pre-reveal visibility', 'RMSE oracle', 'native replay',
               'unsupported inputs', 'tamper rejection'],
    'biologicalSuccessThreshold': 'No biological qualification claimed; retain all gains and losses'})
assay, source, bulk, eligible = load_assay(a.assay)
original = ad.read_h5ad(a.assay.parent / 'original-kang.h5ad')
rows = read(a.assay.parent / 'preparation.json')['originalObservationIndices']
x = original.X[rows].tocsr()
m = bulk['matrix']; counts = np.zeros((len(bulk['groups']), len(bulk['featureIDs'])), dtype=np.uint64)
for i, g in enumerate(bulk['groups']):
    counts[i, m['featureIndices'][m['rowOffsets'][i]:m['rowOffsets'][i+1]]] = m['counts'][m['rowOffsets'][i]:m['rowOffsets'][i+1]]
    np.testing.assert_array_equal(counts[i], np.asarray(x[g['sourceCellIndices']].sum(axis=0)).ravel())
np.testing.assert_array_equal(bulk['featureIDs'], original.var_names)
logs = np.log1p(counts / counts.sum(axis=1)[:,None] * 1e6)
results = []
for specimen in cat['specimens']:
    donor = specimen['donor']
    run = predict(a.assay, a.binary, a.output / 'experiments', donor, assay['observationHours'], assay['intervention'])
    before = summary(run); assert not before['revealed'] and 'comparison' not in before
    train = [d for d in eligible if d != donor]
    c = np.array([eligible[d][0] for d in train]); t = np.array([eligible[d][1] for d in train]); q, target = eligible[donor]
    response = logs[t] - logs[c]; mean = response.mean(axis=0)
    pair_counts = counts[c] + counts[t]
    selected = (pair_counts.sum(axis=0) >= 10) & ((pair_counts > 0).sum(axis=0) >= 2)
    context = logs[c][:,selected]; center = context.mean(axis=0); scale = context.std(axis=0)
    constant = np.all(context == context[0],axis=0); center[constant] = context[0,constant]; scale[constant] = 1
    z = (context-center)/scale/np.sqrt(selected.sum()); query = (logs[q,selected]-center)/scale/np.sqrt(selected.sum())
    ridge = mean + (query @ z.T) @ np.linalg.solve(z @ z.T + np.eye(len(c)), response-mean)
    oracle = {'contextRidge':ridge, 'noChange':np.zeros(len(mean)), 'meanResponse':mean, 'medianResponse':np.median(response,axis=0)}
    for estimate in before['prediction']['estimates']:
        expected = np.maximum(0, logs[q] + oracle[estimate['baseline']])
        np.testing.assert_allclose(estimate['predictedTreated'], expected, atol=1e-9, rtol=1e-9)
    result = reveal(run, a.binary)
    for metric in result['metrics']:
        predicted = np.maximum(0,logs[q]+oracle[metric['model']])
        np.testing.assert_allclose(metric['rmse'], np.sqrt(np.mean((predicted-logs[target])**2)), atol=1e-9, rtol=1e-9)
    verify(run, a.binary)
    results.append({'donor':donor, 'run':run.name, 'metrics':result['metrics'], 'verdict':result['verdict']})
    print(donor + ': ' + result['verdict'], flush=True)
# Fail before making an experiment directory for unsupported user selections.
for donor, hours, intervention in [('unknown',6,'IFN-beta'), ('patient_101',24,'IFN-beta'), ('patient_101',6,'unknown')]:
    try: predict(a.assay, a.binary, a.output/'unsupported', donor, hours, intervention)
    except ValueError: pass
    else: raise AssertionError('Unsupported selection accepted')
assert not (a.output/'unsupported').exists()
# A destructive tamper check uses a separate small copy, never an accepted run.
import shutil
bad = a.output/'tampered'; shutil.copytree(run,bad)
(bad/'prediction/folds/000/prediction.json').write_text('{}')
try: verify(bad,a.binary)
except ValueError: pass
else: raise AssertionError('Changed prediction accepted')
write(a.output/'checks.json', {'passed':True,'donors':len(results),'genesPerDonor':len(mean),
    'nativeReplay':True,'independentNumPyOracle':True,'sourceCountsExact':True,
    'unsupportedInputsRejected':True,'tamperedPredictionRejected':True,'results':results})
