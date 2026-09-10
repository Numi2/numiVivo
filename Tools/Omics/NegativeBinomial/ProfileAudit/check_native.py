#!/usr/bin/env python3
"""Independent GLM refits and stable likelihood checks of measured native profiles."""
import argparse,hashlib,json,warnings
from pathlib import Path
import numpy as np
import statsmodels.api as sm
from scipy.special import gammaln
from scipy.optimize import root

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--refine-reference',action='store_true');a=p.parse_args()
requests=json.loads((a.root/'native-input.json').read_text());actual=json.loads((a.root/'native.json').read_text())
assert [r['id'] for r in requests]==[r['id'] for r in actual]
checks=[];failures=[];messages=[];expected_rejections=[]
limits=dict(meansRelative=2e-5,effect=2e-5,standardError=2e-5,logLikelihood=2e-5,coxReid=2e-5)
for request,native in zip(requests,actual):
 y=np.asarray(request['counts'],dtype=float);x=np.asarray(request['design']);offset=np.asarray(request['offsets']);contrast=np.asarray(request['contrast'])
 steps=[np.arange(int(v),dtype=float) for v in y]
 def likelihood(mu,alpha):
  return float(sum(np.log1p(k*alpha).sum()-gammaln(v+1)+v*np.log(m)-(v+1/alpha)*np.log1p(alpha*m) for k,v,m in zip(steps,y,mu)))
 points=list(native['points'])
 if 'profile' in native:points.append(dict(dispersion=native['profile']['fit']['dispersion'],fit=native['profile']['fit'],selectedProfile=True))
 else:failures.append(dict(id=request['id'],stage='native-profile',error=native.get('profileError')))
 errors={k:0.0 for k in limits};native_identity=0.;reference_objectives=[];point_checks=[]
 for point in points:
  alpha=point['dispersion']
  if not 1e-8<=alpha<=100:
   assert 'fit' not in point and 'domain' in point.get('error','')
   expected_rejections.append(dict(id=request['id'],dispersion=alpha,error=point['error']))
   continue
  if 'fit' not in point:
   failures.append(dict(id=request['id'],dispersion=alpha,stage='native-point',error=point.get('error')));continue
  n=point['fit']
  if not n['converged']:
   failures.append(dict(id=request['id'],dispersion=alpha,stage='native-convergence'));continue
  native_mu=np.array(n['means']);w=native_mu/(1+alpha*native_mu)
  ll=likelihood(native_mu,alpha);r=np.linalg.qr(x*np.sqrt(w[:,None]),mode='r')
  cr=ll-np.log(np.abs(np.diag(r))).sum()
  native_identity=max(native_identity,abs(ll-n['logLikelihood']),abs(cr-n['coxReidLogLikelihood']))
  with warnings.catch_warnings(record=True) as caught:
   warnings.simplefilter('always')
   try:fit=sm.GLM(y,x,offset=offset,family=sm.families.NegativeBinomial(alpha=alpha)).fit(maxiter=1000,tol=1e-12)
   except Exception as error:
    failures.append(dict(id=request['id'],dispersion=alpha,stage='reference-fit',error=str(error)));continue
  messages.extend(dict(id=request['id'],dispersion=alpha,message=str(w.message)) for w in caught)
  beta=fit.params
  if a.refine_reference:
   # Refine IRLS by independently solving the coefficient score equations.
   scale=np.sqrt(np.sum(x*x*(fit.mu/(1+alpha*fit.mu))[:,None],axis=0))
   def score_equations(b):
    mu=np.exp(offset+x@b);return (x.T@((y-mu)/(1+alpha*mu)))/scale
   def jacobian(b):
    mu=np.exp(offset+x@b);curvature=(1+alpha*y)*mu/(1+alpha*mu)**2
    return -(x.T@(curvature[:,None]*x))/scale[:,None]
   solution=root(score_equations,beta,jac=jacobian,method='hybr',options={'xtol':1e-11})
   beta=solution.x
  mu=np.exp(offset+x@beta);weights=mu/(1+alpha*mu);score=np.max(np.abs(x.T@((y-mu)/(1+alpha*mu)))/np.sqrt(np.sum(x*x*weights[:,None],axis=0)))
  if not np.isfinite(score) or score>1e-6:
   failures.append(dict(id=request['id'],dispersion=alpha,stage='reference-score',score=float(score)));continue
  rr=np.linalg.qr(x*np.sqrt(weights[:,None]),mode='r');se=float(np.linalg.norm(np.linalg.solve(rr.T,contrast)))
  ref_ll=likelihood(mu,alpha);ref_cr=float(ref_ll-np.log(np.abs(np.diag(rr))).sum());effect=float(contrast@beta)
  values=dict(meansRelative=float(np.max(np.abs(mu-native_mu)/np.maximum(1,mu))),effect=abs(effect-n['effect']),standardError=abs(se-n['standardError']),logLikelihood=abs(ref_ll-n['logLikelihood']),coxReid=abs(ref_cr-n['coxReidLogLikelihood']))
  for key in errors:errors[key]=max(errors[key],values[key])
  reference_objectives.append(ref_cr);point_checks.append(dict(dispersion=alpha,referenceCoxReid=ref_cr,referenceScaledScore=float(score),errors=values,selectedProfile=point.get('selectedProfile',False)))
  if any(values[k]>limits[k] for k in limits):failures.append(dict(id=request['id'],dispersion=alpha,stage='numerical-comparison',errors=values))
 gap=float(max(reference_objectives)-native['profile']['objective']) if reference_objectives and 'profile' in native else None
 if gap is not None and gap>2e-5:failures.append(dict(id=request['id'],stage='profile-missed-grid-objective',gap=gap))
 checks.append(dict(id=request['id'],points=len(point_checks),maxErrors=errors,stableNativeIdentityError=native_identity,independentGridOverNativeProfile=gap,pointChecks=point_checks))
 (a.root/'native-checks.json').write_text(json.dumps(dict(status='partial',checks=checks,failures=failures,warnings=messages,limits=limits),indent=2,allow_nan=False)+'\n')
 print(json.dumps(dict(id=request['id'],points=len(point_checks),maxCoxReidError=errors['coxReid'],gridGap=gap,failuresSoFar=len(failures))),flush=True)
result=dict(status='passed-selected-valid-profile-numerics' if not failures else 'completed-with-unqualified-points',checks=checks,failures=failures,warnings=messages,expectedDomainRejections=expected_rejections,referenceRefinement=a.refine_reference,limits=limits,
 inputSHA256=hashlib.sha256((a.root/'native-input.json').read_bytes()).hexdigest(),nativeSHA256=hashlib.sha256((a.root/'native.json').read_bytes()).hexdigest(),
 qualification='Selected measured profiles and finite point grid; no global-optimality, calibrated inference or power qualification.')
(a.root/'native-checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(status=result['status'],profiles=len(checks),points=sum(c['points'] for c in checks),failures=len(failures),warnings=len(messages))))
