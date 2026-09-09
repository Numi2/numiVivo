#!/usr/bin/env python3
"""Conditioned NB solver qualification on all eligible real Kang pseudobulks.

This isolates the numerical owner. Fixed dispersion 0.15 is an explicit
validation input, not an estimated biological dispersion or a DE claim.
"""
import argparse,hashlib,json,subprocess,time,warnings
from pathlib import Path
from importlib.metadata import version
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--kang-result',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
counts=pd.read_csv(a.kang_result/'reference-pseudobulk.tsv',sep='\t',index_col=0)
source=json.loads((a.kang_result/'native-report.json').read_text())
design=source['contrasts'][0]['design']
x=np.array(design['rows']);offset=np.log(design['sizeFactorValues']);contrast=np.array(design['contrast'])
assert counts.shape==(16,8894)
assert list(counts.index)==[g['donorID']+'__'+g['condition'] for g in design['observations']]
request=dict(design=x.tolist(),offsets=offset.tolist(),contrast=contrast.tolist(),counts=counts.T.to_numpy().tolist(),dispersion=0.15,profileCount=128)
(a.out/'input.json').write_text(json.dumps(request))
start=time.perf_counter()
r=subprocess.run([str(a.binary),str(a.out/'input.json'),str(a.out/'native.json')],capture_output=True,text=True)
(a.out/'native.log').write_text(r.stdout+r.stderr)
assert r.returncode==0,r.stderr
native=json.loads((a.out/'native.json').read_text())
native_seconds=time.perf_counter()-start
errors=dict(effect=0.,standardError=0.,logLikelihood=0.,meansRelative=0.,coxReid=0.)
reference=[]
diagnostics=[]
rank_deficient=0
for i,gene in enumerate(counts.columns):
    y=counts[gene].to_numpy()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit=sm.GLM(y,x,offset=offset,family=sm.families.NegativeBinomial(alpha=0.15)).fit(maxiter=200,tol=1e-11)
    diagnostics.extend(dict(gene=gene,message=str(w.message),category=w.category.__name__) for w in caught)
    n=native['fits'][i]
    assert n['converged'],(gene,n['maximumScaledScore'])
    deficient=np.linalg.matrix_rank(x[y>0])<x.shape[1]
    assert deficient==n['positiveCountDesignRankDeficient'],gene
    if deficient:
        rank_deficient+=1
        assert 'effect' not in n and 'standardError' not in n and 'coxReidLogLikelihood' not in n
        continue
    effect=float(contrast@fit.params); se=float(np.sqrt(contrast@fit.cov_params()@contrast))
    w=fit.mu/(1+0.15*fit.mu)
    cr=fit.llf-0.5*np.linalg.slogdet(x.T@(w[:,None]*x))[1]
    errors['effect']=max(errors['effect'],abs(effect-n['effect']))
    errors['standardError']=max(errors['standardError'],abs(se-n['standardError']))
    errors['logLikelihood']=max(errors['logLikelihood'],abs(fit.llf-n['logLikelihood']))
    errors['coxReid']=max(errors['coxReid'],abs(cr-n['coxReidLogLikelihood']))
    errors['meansRelative']=max(errors['meansRelative'],float(np.max(np.abs(fit.mu-n['means'])/np.maximum(1,fit.mu))))
    reference.append(dict(gene=gene,effect=effect,standardError=se,logLikelihood=fit.llf))
pd.DataFrame(reference).to_csv(a.out/'statsmodels.tsv',sep='\t',index=False)
# An independent optimizer/refit for the first 128 source-order eligible genes.
# Compare objective value, not arbitrary dispersion location on flat boundaries.
profile_errors=[]
profile_rejections=[]
near_poisson_unqualified=[]
for i,attempt in enumerate(native['profiles']):
    if 'error' in attempt:
        assert np.linalg.matrix_rank(x[counts.iloc[:,i].to_numpy()>0])<x.shape[1]
        assert 'positive-count support is rank deficient' in attempt['error']
        profile_rejections.append(dict(gene=str(counts.columns[i]),reason=attempt['error']))
        continue
    n=attempt['fit']
    y=counts.iloc[:,i].to_numpy()
    def objective(t):
        alpha=np.exp(t)
        f=sm.GLM(y,x,offset=offset,family=sm.families.NegativeBinomial(alpha=alpha)).fit(maxiter=200,tol=1e-11)
        w=f.mu/(1+alpha*f.mu)
        return -(f.llf-0.5*np.linalg.slogdet(x.T@(w[:,None]*x))[1])
    # scipy's gammaln loses precision near alpha=1e-8; evaluate native-selected
    # interior points and the independent optimum only above 1e-5 here.
    opt=minimize_scalar(objective,bounds=(np.log(1e-5),np.log(100)),method='bounded',options={'xatol':1e-7})
    candidates=[(opt.x,opt.fun),(np.log(1e-5),objective(np.log(1e-5))),(np.log(100),objective(np.log(100)))]
    t,negative=min(candidates,key=lambda v:v[1])
    if n['fit']['dispersion']>=1e-5:
        error=abs(-negative-n['objective']);assert error<1e-4,(counts.columns[i],error)
        profile_errors.append(error)
    else:
        near_poisson_unqualified.append(dict(gene=str(counts.columns[i]),dispersion=n['fit']['dispersion'],lowerBoundary=n['lowerBoundary'],reason='Reference log-gamma precision does not qualify this near-Poisson profile'))
thresholds=dict(effect=2e-5,standardError=2e-5,logLikelihood=1e-5,meansRelative=2e-5,coxReid=1e-5)
report=dict(status='passed-conditioned-NB-numerics' if all(errors[k]<thresholds[k] for k in errors) else 'failed',
    genes=counts.shape[1],fullPositiveSupportGenes=len(reference),rankDeficientSupportGenes=rank_deficient,referenceWarnings=diagnostics,profileRejections=profile_rejections,nearPoissonUnqualified=near_poisson_unqualified,observations=16,designColumns=x.shape[1],fixedDispersion=0.15,
    maxErrors=errors,acceptance=thresholds,profileGenes=128,interiorProfileComparisons=len(profile_errors),
    maxInteriorProfileObjectiveError=max(profile_errors,default=None),nativeSeconds=native_seconds,
    versions={k:version(k) for k in ['statsmodels','scipy','numpy','pandas']},
    inputSHA256=hashlib.sha256((a.out/'input.json').read_bytes()).hexdigest(),
    binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
    qualificationBoundary='Conditional GLM numerical qualification on real pseudobulks; no learned dispersion trend, cohort NB DE, Wald calibration or biological model qualification')
(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
assert report['status']=='passed-conditioned-NB-numerics'
