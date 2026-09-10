#!/usr/bin/env python3
"""Independent conditional NB coefficient, information, MAP and Wald checks."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import root
from scipy.stats import nbinom,norm
from statsmodels.stats.multitest import multipletests
TOL=dict(scaledScore=1e-7,effect=2e-6,relativeMean=2e-6,standardError=2e-6,logLikelihood=2e-5,pValue=1e-7)
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();protocol=json.loads((a.root/'protocol.json').read_text());records=[];failures=[];start=time.monotonic()
for case in protocol['cases']:
 d=a.root/case['id'];path=d/'native/report.json.gz'
 if not path.exists():continue
 report=json.loads(gzip.decompress(path.read_bytes()));cohort=report['contrasts'][0];design=cohort['design'];x=np.array(design['rows']);c=np.array(design['contrast']);offsets=np.log(design['sizeFactorValues'])
 counts=pd.read_csv(d/'reference-input/counts.tsv',sep='\t',index_col=0);assert counts.columns.tolist()==[o['sampleIDs'][0] for o in design['observations']]
 rows=[];errors=[]
 for f,diag in zip(cohort['features'],cohort['negativeBinomial']['features'],strict=True):
  if f['status']!='tested':continue
  gene=f['featureID'];y=counts.loc[gene].to_numpy(float);alpha=diag['finalDispersion'];assert 'effectShrinkageError' not in diag
  initial=np.linalg.lstsq(x,np.log(y+.5)-offsets,rcond=None)[0]
  attempts=[]
  for kind in ['MLE','MAP']:
   actual=diag['finalFit'] if kind=='MLE' else diag['effectShrinkageFit'];assert actual['converged']
   precision=0 if kind=='MLE' else 1/np.log(2)**2
   def score(b):
    mu=np.exp(offsets+x@b);return x.T@((y-mu)/(1+alpha*mu))-precision*c*(c@b)
   def info(b):
    mu=np.exp(offsets+x@b);w=(1+alpha*y)*mu/(1+alpha*mu)**2;return x.T@(w[:,None]*x)+precision*np.outer(c,c)
   def score_scale(b):
    if kind=='MAP':return np.sqrt(np.diag(info(b)))
    mu=np.exp(offsets+x@b);return np.sqrt(np.sum(x*x*(mu/(1+alpha*mu))[:,None],axis=0))
   # Both references start from an independent least-squares initialization,
   # never from native coefficients or from a native penalized estimate.
   ref=root(score,initial,jac=lambda b:-info(b),method='hybr',options={'xtol':1e-10});b=ref.x;mu=np.exp(offsets+x@b);scaled=float(np.max(np.abs(score(b))/score_scale(b)))
   actual_beta=np.array(actual['coefficients']);native_score=float(np.max(np.abs(score(actual_beta))/score_scale(actual_beta)))
   effect=float(c@b);information=x.T@((mu/(1+alpha*mu))[:,None]*x) if kind=='MLE' else info(b);sd=float(np.sqrt(c@np.linalg.solve(information,c)));actual_sd=actual['standardError'] if kind=='MLE' else actual['posteriorStandardDeviation']
   metrics=dict(scaledScore=max(scaled,native_score),effect=abs(effect-actual['effect']),relativeMean=float(np.max(np.abs(mu-np.array(actual['means']))/mu)),standardError=abs(sd-actual_sd),logLikelihood=abs(float(nbinom.logpmf(y,1/alpha,1/(1+alpha*mu)).sum())-actual['logLikelihood']),pValue=abs(float(2*norm.sf(abs(effect/sd)))-f['pValue']) if kind=='MLE' else 0.)
   failed=[key for key,v in metrics.items() if not np.isfinite(v) or v>TOL[key]]
   if kind=='MLE':
    if abs(f['log2FoldChange']-effect/np.log(2))>2e-6:failed.append('log2Effect')
    if abs(f['standardError']-sd/np.log(2))>2e-6:failed.append('log2SE')
   attempt=dict(kind=kind,metrics=metrics,referenceSuccess=bool(ref.success),referenceMessage=str(ref.message),failed=failed);attempts.append(attempt)
   if failed:errors.append(dict(featureID=gene,kind=kind,failed=failed,metrics=metrics))
  rows.append(dict(featureIndex=f['featureIndex'],featureID=gene,attempts=attempts))
 tested=[f for f in cohort['features'] if f['status']=='tested'];np.testing.assert_allclose(multipletests([f['pValue'] for f in tested],method='fdr_bh')[1],[f['adjustedPValue'] for f in tested],rtol=2e-13,atol=1e-15)
 summary=dict(case=case['id'],cellGroup=case['cellGroup'],testedGenes=len(rows),conditionalFits=2*len(rows),failures=errors,referenceUnsuccessfulFlags=sum(not t['referenceSuccess'] for r in rows for t in r['attempts']),maximumErrors={k:max(t['metrics'][k] for r in rows for t in r['attempts']) for k in TOL},independentBHVerified=True)
 records.append(summary);failures.extend(errors);(d/'model-checks.json.gz').write_bytes(gzip.compress(json.dumps(rows,sort_keys=True,allow_nan=False).encode(),mtime=0));(d/'model-check.json').write_text(json.dumps(summary,sort_keys=True,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='failures'}),flush=True)
result=dict(status='passed-conditional-model-checks' if not failures else 'failed',tolerances=TOL,failedFits=len(failures),cases=records,seconds=time.monotonic()-start,qualification='Fixed final dispersion coefficient/information/MAP and Wald/BH arithmetic checks, not dispersion calibration, posterior coverage or FDR qualification',checkerSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(a.root/'model-check.json').write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n');assert not failures
