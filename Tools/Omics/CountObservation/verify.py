"""Independent adaptive quadrature reference, including high precision boundaries."""
from pathlib import Path
import json,sys,math,time,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import gamma
import mpmath as mp
cases=json.loads(Path(sys.argv[1]).read_bytes());records=[json.loads(x) for x in Path(sys.argv[2]).read_text().splitlines()]
out=Path(sys.argv[3]);assert not out.exists()
assert len(cases)==len(records)
checks=[];failures=[];started=time.time()
for index,(c,native) in enumerate(zip(cases,records)):
 assert c['id']==native['id']
 if 'error' in native:
  failures.append({'id':c['id'],'nativeError':native['error']});continue
 y=np.array(c['counts'],float);e=np.array(c['libraryCounts'],float)/1e6
 a,b,phi=c['gammaPriorShape'],c['gammaPriorRatePerCPM'],c['cellDispersion'];shape=float(y.sum())+a
 # Independent score-root construction and QUADPACK adaptive integrals. Native
 # uses centered Taylor arithmetic and composite 16-point Gauss-Legendre panels.
 def score(x):
  if phi==0:return shape-(b+e.sum())*np.exp(x)
  q=expit(np.log(phi*e)+x)
  return a-b*np.exp(x)+float(np.sum(y*(1-q)-q/phi))
 mode=brentq(score,-690,690,xtol=1e-13)
 if c['id'].startswith('boundary'):
  mp.mp.dps=60
  yy=[mp.mpf(int(v)) for v in c['counts']];ee=[mp.mpf(int(v))/10**6 for v in c['libraryCounts']]
  aa,bb,pp=map(lambda v:mp.mpf(str(v)),[a,b,phi]);mm=mp.mpf(mode)
  def logpost(x):
   if pp==0:return (sum(yy)+aa)*x-(bb+sum(ee))*mp.exp(x)
   return (sum(yy)+aa)*x-bb*mp.exp(x)-sum((v+1/pp)*mp.log1p(pp*l*mp.exp(x)) for v,l in zip(yy,ee))
  baseline=logpost(mm)
  def density(d):return float(mp.exp(logpost(mm+mp.mpf(d))-baseline))
 else:
  def logpost(x):
   if phi==0:return shape*x-(b+e.sum())*np.exp(x)
   return shape*x-b*np.exp(x)-float(np.sum((y+1/phi)*np.logaddexp(0,np.log(phi*e)+x)))
  baseline=logpost(mode)
  def density(d):return math.exp(logpost(mode+d)-baseline)
 q=expit(np.log(phi*e)+mode) if phi else None
 curvature=(b+e.sum())*math.exp(mode) if phi==0 else b*math.exp(mode)+float(np.sum((y+1/phi)*q*(1-q)))
 width=1/math.sqrt(curvature)
 # Work in standard-width coordinates so concentrated posteriors are not missed.
 left=-1.;right=1.
 while density(left*width)>1e-28:left*=1.4
 while density(right*width)*math.exp(min(600,2*right*width))>1e-28:right*=1.4
 def integrate(f,l=left,h=right):
  # Split at the mode; adaptive integration is independent of native partitions.
  intervals=[(l,0),(0,h)] if l<0<h else [(l,h)]
  return sum(quad(lambda z:density(z*width)*f(z*width),lo,hi,epsabs=2e-11,epsrel=2e-10,limit=300)[0] for lo,hi in intervals)
 z=integrate(lambda d:1)
 rate=math.exp(mode);logcenter=float(np.logaddexp(0,mode))
 dr=integrate(math.expm1)/z
 vr=integrate(lambda d:(math.expm1(d)-dr)**2)/z
 def dl(d):return math.log1p(float(expit(mode))*math.expm1(d)) if abs(d)<0.5 else float(np.logaddexp(0,mode+d))-logcenter
 ml=integrate(dl)/z;vl=integrate(lambda d:(dl(d)-ml)**2)/z
 quantiles=[]
 for p in [.025,.975]:
  root=brentq(lambda x:integrate(lambda d:1,left,x)/z-p,left,right,xtol=1e-11)
  quantiles.append(float(np.logaddexp(0,mode+root*width)))
 expected={'meanCPM':rate*(1+dr),'varianceCPM':rate*rate*vr,'meanLog1pCPM':logcenter+ml,'varianceLog1pCPM':vl,'lowerLog1pCPM':quantiles[0],'upperLog1pCPM':quantiles[1]}
 p=native['posterior'];errors={k:abs(p[k]-v)/max(1e-12,abs(v)) for k,v in expected.items()}
 passed=all(err<2e-7 for err in errors.values()) and all(abs(p[k]-expected[k])<2e-7 for k in ['lowerLog1pCPM','upperLog1pCPM'])
 # Analytic Gamma-Poisson posterior is an additional, quadrature-free check.
 analytic=None
 if phi==0:
  mean=shape/(b+e.sum());variance=shape/(b+e.sum())**2
  bounds=np.log1p(gamma.ppf([.025,.975],shape,scale=1/(b+e.sum())))
  analytic=max(abs(p['meanCPM']/mean-1),abs(p['varianceCPM']/variance-1),max(abs(np.array([p['lowerLog1pCPM'],p['upperLog1pCPM']])-bounds)))
  passed &= analytic<2e-7
 depths=np.array(c['plannedLibraryCounts'])/1e6;E=depths.sum();pred=native['prediction']
 future={'meanGeneCounts':E*expected['meanCPM'],'conditionalPoissonVariance':E*expected['meanCPM'],'conditionalCellOverdispersionVariance':phi*float(depths@depths)*(expected['varianceCPM']+expected['meanCPM']**2),'latentRateVariance':E*E*expected['varianceCPM']}
 future['totalGeneCountVariance']=future['conditionalPoissonVariance']+future['conditionalCellOverdispersionVariance']+future['latentRateVariance']
 ferr={k:abs(pred[k]-v)/max(1e-12,abs(v)) for k,v in future.items()};passed &= max(ferr.values())<2e-7
 entry={'id':c['id'],'passed':bool(passed),'maximumRelativePosteriorError':max(errors.values()),'maximumRelativePredictiveMomentError':max(ferr.values()),'analyticGammaError':analytic,'reference':expected}
 checks.append(entry)
 if not passed:failures.append(entry)
 if (index+1)%100==0:print(json.dumps({'checked':index+1,'failed':len(failures)}),flush=True)
result={'status':'passed' if not failures else 'failed','cases':len(cases),'checked':len(checks),'failures':failures,'checks':checks,'seconds':time.time()-started,'inputSHA256':hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest(),'nativeSHA256':hashlib.sha256(Path(sys.argv[2]).read_bytes()).hexdigest(),'referenceRecipeSHA256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['checks','failures']}));sys.exit(bool(failures))
