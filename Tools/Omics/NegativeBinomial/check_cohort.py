#!/usr/bin/env python3
"""Independent checks of native cohort trend, MAP and Wald calculations.

Consumes the full real Kang CLI benchmark, not synthetic scientific evidence.
Statsmodels supplies coefficient fits; SciPy supplies a separate MAP optimizer.
"""
import argparse,hashlib,json,warnings
from pathlib import Path
from importlib.metadata import version
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import minimize,minimize_scalar
from scipy.special import polygamma,erfc
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--kang-result',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
report=json.loads((a.kang_result/'native-report.json').read_text())
c=report['contrasts'][0];d=c['design'];nb=c['negativeBinomial'];trend=nb['trend']
x=np.asarray(d['rows']);offset=np.log(d['sizeFactorValues']);contrast=np.asarray(d['contrast'])
counts=pd.read_csv(a.kang_result/'reference-pseudobulk.tsv',sep='\t',index_col=0)
assert list(counts.index)==[o['donorID']+'__'+o['condition'] for o in d['observations']]
features=c['features'];diagnostics=nb['features'];indices=trend['referenceFeatureIndices']
means=np.array([features[i]['meanNormalizedCount'] for i in indices])
alpha=np.array([diagnostics[i]['geneWiseDispersion'] for i in indices])
target=trend['intercept']+trend['inverseMeanCoefficient']/means
residual=np.log(alpha/target);center=np.median(residual)
robust=float((np.median(abs(residual-center))/0.6744897501960817)**2)
sampling=float(polygamma(1,d['residualDegreesOfFreedom']/2))
prior=max(c['request']['negativeBinomialOptions'].get('minimumPriorVariance',0.25),robust-sampling)
assert abs(robust-trend['robustLogResidualVariance'])<1e-10
assert abs(sampling-trend['samplingLogVariance'])<1e-10
assert abs(prior-trend['priorLogVariance'])<1e-10
trend_gap=None
if trend['method']=='parametric':
    scale=np.median(means);inv=scale/means
    def loss(theta):
        t=theta[0]+theta[1]*inv
        if np.any(t<=0):return np.inf
        r=abs(np.log(alpha/t))
        return float(np.sum(np.where(r<=1.345,r*r/2,1.345*(r-1.345/2))))
    initial=np.array([np.exp(np.median(np.log(alpha))),np.exp(np.median(np.log(alpha)))*.1])
    opt=minimize(loss,initial,method='Nelder-Mead',bounds=[(0,None),(0,None)],options={'maxiter':3000,'xatol':1e-10,'fatol':1e-10})
    assert opt.success,opt.message
    trend_gap=loss([trend['intercept'],trend['inverseMeanCoefficient']/scale])-opt.fun
    assert abs(trend_gap)<1e-5,trend_gap
else:assert abs(trend['intercept']-np.exp(np.median(np.log(alpha))))<1e-10
errors=dict(effect=0.,standardError=0.,z=0.,p=0.)
profile_errors=[];reference_warnings=[];tested=0
for i,f in enumerate(features):
    if f['status']!='tested':
        assert 'pValue' not in f and 'adjustedPValue' not in f
        continue
    tested+=1
    y=counts[f['featureID']].to_numpy();a_final=diagnostics[i]['finalDispersion']
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit=sm.GLM(y,x,offset=offset,family=sm.families.NegativeBinomial(alpha=a_final)).fit(maxiter=200,tol=1e-11)
    reference_warnings.extend(dict(gene=f['featureID'],message=str(w.message)) for w in caught)
    effect=float(contrast@fit.params);se=float(np.sqrt(contrast@fit.cov_params()@contrast));z=effect/se
    errors['effect']=max(errors['effect'],abs(effect/np.log(2)-f['log2FoldChange']))
    errors['standardError']=max(errors['standardError'],abs(se/np.log(2)-f['standardError']))
    errors['z']=max(errors['z'],abs(z-f['zStatistic']))
    errors['p']=max(errors['p'],abs(float(erfc(abs(z)/np.sqrt(2)))-f['pValue']))
    if len(profile_errors)<32 and not diagnostics[i]['dispersionOutlier']:
        target=diagnostics[i]['trendDispersion']
        def objective(t):
            alpha=np.exp(t)
            ref=sm.GLM(y,x,offset=offset,family=sm.families.NegativeBinomial(alpha=alpha)).fit(maxiter=200,tol=1e-11)
            w=ref.mu/(1+alpha*ref.mu)
            return -(ref.llf-.5*np.linalg.slogdet(x.T@(w[:,None]*x))[1]-(t-np.log(target))**2/(2*prior))
        opt=minimize_scalar(objective,bounds=(np.log(1e-5),np.log(100)),method='bounded',options={'xatol':1e-7})
        assert opt.success
        error=abs(objective(np.log(a_final))-opt.fun)
        assert error<1e-4,(f['featureID'],error)
        profile_errors.append(dict(gene=f['featureID'],objectiveError=error))
assert tested>0 and len(profile_errors)==32
assert all(v<2e-5 for v in errors.values()),errors
# Independent BH on exactly the declared available hypotheses.
from statsmodels.stats.multitest import multipletests
rows=[f for f in features if f['status']=='tested']
q=multipletests([f['pValue'] for f in rows],method='fdr_bh')[1]
assert np.max(abs(q-[f['adjustedPValue'] for f in rows]))<1e-12
summary=dict(status='passed-native-cohort-calculations',testedGenes=tested,trendMethod=trend['method'],
    trendReferenceGenes=len(indices),trendObjectiveGap=trend_gap,priorLogVariance=prior,maxWaldErrors=errors,
    mapChecks=profile_errors,referenceWarnings=reference_warnings,exactBH=True,
    sourceReportSHA256=hashlib.sha256((a.kang_result/'native-report.json').read_bytes()).hexdigest(),
    versions={k:version(k) for k in ['statsmodels','scipy','numpy','pandas']},
    boundary='Numerical checks conditional on the native model and one real study; not independent biological validity or repeated-sampling calibration')
(a.out/'report.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
