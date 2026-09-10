#!/usr/bin/env python3
"""Independent R empirical-prior and SciPy conditional score/Laplace checks."""
import argparse,gzip,hashlib,json,platform,time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import root
from scipy.stats import nbinom
from scipy.sparse import csr_matrix

TOL=dict(scaledScore=1e-7,effect=2e-6,relativeMean=2e-6,posteriorSD=2e-6,logLikelihood=2e-5,objective=2e-5)
def load(p): return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--case');a=p.parse_args()
start=time.monotonic();protocol=load(a.root/'protocol.json');studies=[];all_failures=[]
prior_references={r['case']:r for r in load(a.root/'prior-reference.json')['records']}
for rec in protocol['records']:
 if a.case is not None and rec['id']!=a.case:continue
 d=a.root/rec['id'];old=load(Path(rec['baselineReport']));new=load(d/'native/report.json.gz')
 assert sha(Path(rec['baselineReport']))==rec['baselineReportSHA256']
 ob=old['contrasts'][0];nb=new['contrasts'][0];errors=[];rows=[]
 for field in ['metadata','pseudobulk','quality','canonicalNonzeros']:assert old[field]==new[field],field
 for field in ['features','testedFeatures','design','method','multiplicityScope']:assert ob[field]==nb[field],field
 assert ob['negativeBinomial']['trend']==nb['negativeBinomial']['trend']
 old_request=dict(nb['request']);old_request['negativeBinomialOptions']=dict(old_request['negativeBinomialOptions']);del old_request['negativeBinomialOptions']['effectPriorEstimation'];old_request['negativeBinomialOptions']['effectPriorStandardDeviationLog2']=1;assert old_request==ob['request']
 estimate=nb['negativeBinomial']['effectPriorEstimate'];ref_prior=prior_references[rec['id']];assert 'effectPriorEstimationError' not in nb['negativeBinomial']
 assert estimate['referenceFeatureIndices']==ref_prior['referenceFeatureIndices']
 eligible=[f['featureIndex'] for f,dg in zip(nb['features'],nb['negativeBinomial']['features'],strict=True) if f['status']=='tested' and 'supportResolution' not in dg]
 assert estimate['eligibleFeatureIndices']==eligible
 reference_set=set(ref_prior['referenceFeatureIndices'])
 assert estimate['excludedHighEffectFeatureIndices']==[i for i in eligible if i not in reference_set]
 assert estimate['method']=='weightedUpperQuantile'
 assert abs(estimate['absoluteEffectQuantileLog2']/1.959963984540054-estimate['unboundedStandardDeviationLog2'])<1e-14
 assert estimate['standardDeviationFloorReached']==(ref_prior['rawSDLog2']<=.01)
 for native_key,r_key in [('unboundedStandardDeviationLog2','rawSDLog2'),('priorStandardDeviationLog2','priorSDLog2'),('weightSum','weightSum'),('effectiveReferenceFeatures','effectiveReferenceFeatures')]:
  assert abs(estimate[native_key]-ref_prior[r_key])<=2e-12*max(1,abs(ref_prior[r_key])),(rec['id'],native_key,estimate[native_key],ref_prior[r_key])
 m=new['pseudobulk']['matrix'];counts=csr_matrix((m['counts'],m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount'])).toarray().T
 design=nb['design'];x_full=np.array(design['rows']);c_full=np.array(design['contrast']);offsets_full=np.log(design['sizeFactorValues']);source_rows=design['sourcePseudobulkIndices']
 for before,diag,feature in zip(ob['negativeBinomial']['features'],nb['negativeBinomial']['features'],nb['features'],strict=True):
  original={k:v for k,v in diag.items() if k not in ['effectShrinkageFit','effectShrinkageError']};assert original=={k:v for k,v in before.items() if k not in ['effectShrinkageFit','effectShrinkageError']}
  if feature['status']!='tested':assert 'effectShrinkageFit' not in diag and 'effectShrinkageError' not in diag;continue
  gene=feature['featureIndex'];fit=diag.get('effectShrinkageFit');assert fit and fit['converged'] and 'effectShrinkageError' not in diag
  s=diag.get('supportResolution');retained=s['retainedObservationIndices'] if s else list(range(len(x_full)))
  x=np.array(s['rows']) if s else x_full;c=np.array(s['contrast']) if s else c_full
  offsets=offsets_full[retained];y=counts[gene,source_rows][retained].astype(float);alpha=diag['finalDispersion'];prior=np.log(2)*estimate['priorStandardDeviationLog2'];precision=1/prior**2
  assert fit['dispersion']==alpha and abs(fit['priorStandardDeviation']-prior)<1e-14
  def score(beta):
   mu=np.exp(offsets+x@beta)
   return x.T@((y-mu)/(1+alpha*mu))-c*(c@beta)*precision
  def info(beta):
   mu=np.exp(offsets+x@beta);w=(1+alpha*y)*mu/(1+alpha*mu)**2
   return x.T@(w[:,None]*x)+np.outer(c,c)*precision
  initial=np.array(diag['finalFit']['coefficients'])
  reference=root(score,initial,jac=lambda b:-info(b),method='hybr',options={'xtol':1e-10})
  beta=reference.x;mu=np.exp(offsets+x@beta);effect=float(c@beta);sd=float(np.sqrt(c@np.linalg.solve(info(beta),c)))
  scaled=float(np.max(np.abs(score(beta))/np.sqrt(np.diag(info(beta)))))
  ll=float(nbinom.logpmf(y,1/alpha,1/(1+alpha*mu)).sum());obj=ll-effect**2*precision/2
  actual=np.array(fit['coefficients']);actual_info=info(actual);actual_score=float(np.max(np.abs(score(actual))/np.sqrt(np.diag(actual_info))))
  metrics=dict(scaledScore=max(scaled,actual_score,fit['maximumScaledScore']),effect=abs(effect-fit['effect']),relativeMean=float(np.max(np.abs(mu-np.array(fit['means']))/mu)),posteriorSD=abs(sd-fit['posteriorStandardDeviation']),logLikelihood=abs(ll-fit['logLikelihood']),objective=abs(obj-fit['objective']))
  failed=[key for key,value in metrics.items() if not np.isfinite(value) or value>TOL[key]]
  initial_effect=float(c@initial);mle_objective=diag['finalFit']['logLikelihood']-initial_effect**2*precision/2
  if fit['objective']<mle_objective-TOL['objective']:failed.append('penalizedObjectiveImprovement')
  if abs(fit['effect'])>abs(initial_effect)+TOL['effect'] or fit['effect']*initial_effect < -TOL['effect']:failed.append('effectTowardZero')
  if fit['posteriorStandardDeviation']>prior+TOL['posteriorSD']:failed.append('posteriorSDPriorBound')
  if failed:errors.append(dict(featureIndex=gene,featureID=feature['featureID'],failed=failed,metrics=metrics,referenceMessage=str(reference.message)))
  rows.append(dict(featureIndex=gene,featureID=feature['featureID'],retainedObservations=len(y),effectLog2=fit['effect']/np.log(2),posteriorSDLog2=fit['posteriorStandardDeviation']/np.log(2),unshrunkEffectLog2=feature['log2FoldChange'],referenceSuccess=bool(reference.success),referenceMessage=str(reference.message),metrics=metrics,failed=failed))
 summary=nb['negativeBinomial']['effectShrinkage'];assert summary['eligibleFeatures']==len(rows)==summary['convergedFeatures'] and summary['failedFeatures']==0
 assert summary['priorStandardDeviationLog2']==estimate['priorStandardDeviationLog2']
 record=dict(case=rec['id'],priorEstimate=estimate,RReferencePrior=ref_prior,testedFeatures=len(rows),BH005Calls=sum(f.get('adjustedPValue',1)<=.05 for f in nb['features']),baselineInferenceExactlyEqual=True,baselineCountsDesignAndDiagnosticsExactlyEqual=True,activeSupportFits=sum('supportResolution' in f and 'effectShrinkageFit' in f for f in nb['negativeBinomial']['features']),referenceSolverFlagFalse=sum(not r['referenceSuccess'] for r in rows),maximumErrors={key:max(r['metrics'][key] for r in rows) for key in TOL},failures=errors,reportSHA256=sha(d/'native/report.json.gz'))
 (d/'reference-checks.json.gz').write_bytes(gzip.compress(json.dumps(rows,sort_keys=True,allow_nan=False).encode(),mtime=0));write(d/'check.json',record);studies.append(record);all_failures.extend(errors);print(json.dumps({k:v for k,v in record.items() if k!='failures'}),flush=True)
assert studies,'No selected case found'
write(a.root/('check.json' if a.case is None else 'check-'+a.case+'.json'),dict(status='passed-conditional-numerical-checks' if not all_failures else 'failed',tolerances=TOL,studies=studies,failedFeatures=len(all_failures),seconds=time.monotonic()-start,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,qualification=protocol['qualification'],checkerSHA256=sha(Path(__file__))))
assert not all_failures
