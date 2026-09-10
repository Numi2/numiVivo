#!/usr/bin/env python3
"""Separate stored-parameter arithmetic from independent-refit stopping error.

This diagnostic preserves, and does not override, check_models.py's strict gate.
"""
import argparse, gzip, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import root
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
a = p.parse_args()
original = json.loads((a.root/'model-check.json').read_text())
results = []
for case in original['cases']:
    d = a.root/case['case']
    cohort = json.loads(gzip.decompress((d/'native/report.json.gz').read_bytes()))['contrasts'][0]
    design = cohort['design']
    x = np.array(design['rows']); c = np.array(design['contrast'])
    offsets = np.log(design['sizeFactorValues'])
    counts = pd.read_csv(d/'reference-input/counts.tsv', sep='\t', index_col=0)
    records = []; exceptions = []
    for f, diag in zip(cohort['features'], cohort['negativeBinomial']['features'], strict=True):
        if f['status'] != 'tested': continue
        y = counts.loc[f['featureID']].to_numpy(float); alpha = diag['finalDispersion']
        fit = diag['finalFit']; b0 = np.array(fit['coefficients'])
        def score(b):
            mu = np.exp(offsets+x@b)
            return x.T@((y-mu)/(1+alpha*mu))
        def observed(b):
            mu = np.exp(offsets+x@b)
            return x.T@(((1+alpha*y)*mu/(1+alpha*mu)**2)[:, None]*x)
        def expected(b):
            mu = np.exp(offsets+x@b)
            return x.T@((mu/(1+alpha*mu))[:, None]*x)
        initial = np.linalg.lstsq(x, np.log(y+.5)-offsets, rcond=None)[0]
        ref = root(score, initial, jac=lambda b: -observed(b), method='hybr', options={'xtol': 1e-10})
        scaled_score = float(np.max(np.abs(score(ref.x))/np.sqrt(np.diag(expected(ref.x)))))
        # HYBR's progress flag is distinct from the score equation residual.
        # Preserve every flag, and require stationarity independently.
        assert scaled_score < 1e-9, (f['featureID'], ref.message, scaled_score)
        se = np.sqrt(c@np.linalg.solve(expected(ref.x), c))
        ref_p = float(2*norm.sf(abs(c@ref.x/se)))
        stored_p = float(2*norm.sf(abs(fit['effect']/fit['standardError'])))
        records.append(dict(featureID=f['featureID'], nativeP=f['pValue'], refitP=ref_p,
                            referenceSuccess=bool(ref.success), referenceMessage=str(ref.message),
                            referenceScaledScore=scaled_score,
                            storedArithmeticError=abs(stored_p-f['pValue'])))
        if abs(ref_p-f['pValue']) > original['tolerances']['pValue']:
            mu = np.exp(offsets+x@b0); q = np.linalg.solve(expected(b0), c)
            variance = c@q; effect = c@b0
            variance_gradient = -x.T@((mu/(1+alpha*mu)**2)*(x@q)**2)
            z = effect/np.sqrt(variance)
            z_gradient = c/np.sqrt(variance)-effect*variance_gradient/(2*variance**1.5)
            p_gradient = -2*norm.pdf(abs(z))*np.sign(z)*z_gradient
            predicted = float(p_gradient@np.linalg.solve(observed(b0), score(b0)))
            exceptions.append(dict(featureID=f['featureID'], nativeP=f['pValue'], refitP=ref_p,
                observedChange=ref_p-stored_p, firstOrderNewtonChange=predicted,
                firstOrderResidual=abs(ref_p-stored_p-predicted)))
    native_q = multipletests([r['nativeP'] for r in records], method='fdr_bh')[1]
    ref_q = multipletests([r['refitP'] for r in records], method='fdr_bh')[1]
    summary = dict(case=case['case'], cellGroup=case['cellGroup'], testedGenes=len(records),
        maximumStoredArithmeticError=max(r['storedArithmeticError'] for r in records),
        maximumRefitPValueDifference=max(abs(r['nativeP']-r['refitP']) for r in records),
        maximumRefitBHDifference=float(np.max(np.abs(native_q-ref_q))),
        changedBH005Decisions=int(np.count_nonzero((native_q<.05)!=(ref_q<.05))),
        referenceUnsuccessfulFlags=[r for r in records if not r['referenceSuccess']],
        strictExceptions=exceptions)
    (d/'probability-records.json.gz').write_bytes(gzip.compress(json.dumps(records, sort_keys=True, allow_nan=False).encode(), mtime=0))
    results.append(summary); print(json.dumps(summary), flush=True)
result = dict(status='diagnostic-only-original-strict-gate-retained', cases=results,
    originalStrictGateStatus=original['status'], originalFailedFits=original['failedFits'],
    originalCheckSHA256=hashlib.sha256((a.root/'model-check.json').read_bytes()).hexdigest(),
    diagnosticSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    qualification='Stored-parameter Wald arithmetic and conditional optimization sensitivity; no FDR calibration claim.')
(a.root/'probability-diagnostic.json').write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
assert max(c['maximumStoredArithmeticError'] for c in results) < 1e-14
assert sum(c['changedBH005Decisions'] for c in results) == 0
