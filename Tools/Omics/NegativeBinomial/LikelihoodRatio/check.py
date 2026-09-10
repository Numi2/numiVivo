#!/usr/bin/env python3
"""Independent count-model/null-score/LR/BH checks; prepare pinned edgeR inputs."""
import argparse,gzip,hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.stats import chi2
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--case',help='Check one completed case; full qualification still requires every protocol case');p.add_argument('--family',choices=['null','treatment']);p.add_argument('--refine-reference',action='store_true',help='Recompute cancellation-sensitive reference differences with 80 digits, retaining original discrepancies');a=p.parse_args()
protocol=json.loads((a.root/'protocol.json').read_text())
cases=[r for r in protocol['cases'] if (a.case is None or r['id']==a.case) and (a.family is None or r['family']==a.family)];assert cases
def sha(b):return hashlib.sha256(b).hexdigest()
def gzwrite(p,x):p.write_bytes(gzip.compress(json.dumps(x,allow_nan=False,separators=(',',':')).encode(),mtime=0))
def bh(p):
 p=np.asarray(p);order=np.argsort(p,kind='stable');out=np.empty(len(p));out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return out
for rec in cases:
 d=a.root/rec['id'];run=json.loads((d/'run.json').read_text());raw=(d/'native.json.gz').read_bytes()
 assert sha(raw)==run['outputSHA256'];assert sha(gzip.decompress(raw))==run['logicalSHA256'];native=json.loads(gzip.decompress(raw))
 raw=subprocess.check_output(['ssh',protocol['host'],'cat',rec['source']]);assert sha(raw)==rec['sourceSHA256'];expanded=gzip.decompress(raw);assert sha(expanded)==rec['logicalSHA256'];source=json.loads(expanded)
 old=source['contrasts'][0];design=old['design'];counts=np.zeros((len(source['pseudobulk']['groups']),len(source['metadata']['features'])),dtype=np.uint64);csr=source['pseudobulk']['matrix']
 for i in range(counts.shape[0]):
  start,end=csr['rowOffsets'][i:i+2];counts[i,csr['featureIndices'][start:end]]=csr['counts'][start:end]
 y=counts[design['sourcePseudobulkIndices']];x=np.asarray(design['rows']);offsets=np.log(design['sizeFactorValues']);contrast=np.asarray(design['contrast'])
 for name in ['originalFeaturesExactlyEqual','originalDesignExactlyEqual','originalDiagnosticsExactlyEqual','lrtDesignExactlyEqual','lrtTrendExactlyEqual']:assert native[name],(rec['id'],name)
 request=dict(native['request']);options=dict(request['negativeBinomialOptions']);assert options.pop('testMethod')=='likelihoodRatio';request['negativeBinomialOptions']=options;assert request==old['request']
 metrics={k:0.0 for k in ['scaledScore','relativeMean','constraint','statistic','pValue','BH']};failures=[];unavailable=[];tested=[];r_groups={};details=[];refinements=[]
 for gene,oldfeature,diag in zip(native['genes'],old['features'],old['negativeBinomial']['features'],strict=True):
  f=gene['feature'];idx=f['featureIndex'];assert f['featureID']==oldfeature['featureID']==source['metadata']['features'][idx]['id'];assert gene['fullFitExactlyEqual']
  assert f['totalCounts']==int(y[:,idx].sum()) and f['expressingPseudobulks']==int((y[:,idx]>0).sum())
  if oldfeature['status']!='tested':assert f['status']==oldfeature['status'] and gene.get('likelihoodRatio') is None;continue
  if f['status']!='tested':unavailable.append(dict(featureIndex=idx,error=gene.get('error'),diagnostics=gene.get('likelihoodRatio')));continue
  tested.append(f);fit=gene['likelihoodRatio'];assert fit.get('error') is None
  assert fit['nullConverged'] and fit['degreesOfFreedom']==1
  resolution=diag.get('supportResolution');xx=np.asarray(resolution['rows']) if resolution else x;cc=np.asarray(resolution['contrast']) if resolution else contrast
  rows=resolution['retainedObservationIndices'] if resolution else list(range(len(x)));yy=y[rows,idx].astype(float);oo=offsets[rows];alpha=diag['finalDispersion']
  beta=np.asarray(fit['nullCoefficients']);mu=np.asarray(fit['nullMeans']);full=np.asarray(diag['finalFit']['means'])
  pivot=max(range(len(cc)),key=lambda i:(abs(cc[i]),i));free=[i for i in range(len(cc)) if i!=pivot]
  basis=np.eye(len(cc))[:,free];basis[pivot,:]=-cc[free]/cc[pivot];reduced=xx@basis
  score=reduced.T@((yy-mu)/(1+alpha*mu));info=(reduced**2).T@(mu/(1+alpha*mu));scaled=float(np.max(np.abs(score)/np.sqrt(info)))
  mm=np.exp(xx@beta+oo);constraint=float(abs(cc@beta)/max(1,np.linalg.norm(cc)*np.linalg.norm(beta)))
  mf=full.astype(np.longdouble);mn=mu.astype(np.longdouble);aa=np.longdouble(alpha);yl=yy.astype(np.longdouble)
  statistic=float(2*np.sum(yl*(np.log(mf)-np.log(mn))-(yl+1/aa)*(np.log1p(aa*mf)-np.log1p(aa*mn))))
  probability=float(chi2.sf(max(0,statistic),1));values=dict(scaledScore=scaled,relativeMean=float(np.max(np.abs(mm-mu)/np.maximum(1,mu))),constraint=constraint,statistic=abs(statistic-fit['rawStatistic']),pValue=abs(probability-f['pValue']))
  if a.refine_reference and (values['statistic']>2e-7 or values['pValue']>2e-8):
   import mpmath as mp
   before=dict(values)
   with mp.workdps(80):
    ap=mp.mpf(float(alpha));terms=[]
    for yv,mfv,mnv in zip(yy,full,mu):
     yp=mp.mpf(int(yv));fp=mp.mpf(float(mfv));np_=mp.mpf(float(mnv))
     terms.append(yp*mp.log(fp/np_)-(yp+1/ap)*mp.log((1+ap*fp)/(1+ap*np_)))
    precise=2*mp.fsum(terms);statistic=float(precise);probability=float(mp.erfc(mp.sqrt(max(0,precise)/2)))
   values['statistic']=abs(statistic-fit['rawStatistic']);values['pValue']=abs(probability-f['pValue'])
   refinements.append(dict(featureIndex=idx,initialErrors=before,refinedErrors=dict(values),precisionDigits=80))
  for key,value in values.items():metrics[key]=max(metrics[key],value)
  limits=dict(scaledScore=1.1e-7,relativeMean=1e-10,constraint=1e-10,statistic=2e-6,pValue=2e-7)
  errors=[key for key,v in values.items() if not np.isfinite(v) or v>limits[key]]
  assert f['pValue']==fit['pValue'];assert f['log2FoldChange']==oldfeature['log2FoldChange'] and f['standardError']==oldfeature['standardError'] and f['zStatistic']==oldfeature['zStatistic']
  if errors:failures.append(dict(featureIndex=idx,metrics=values,failed=errors))
  details.append(dict(featureIndex=idx,nativeStatistic=fit['statistic'],referenceStatistic=statistic,referencePValue=probability,**values))
  if rec['family']=='treatment':
   key=json.dumps([xx.tolist(),cc.tolist(),oo.tolist()]);group=r_groups.setdefault(key,dict(design=xx.tolist(),contrast=cc.tolist(),offsets=oo.tolist(),featureIndices=[],counts=[],dispersions=[]))
   group['featureIndices'].append(idx);group['counts'].append(yy.tolist());group['dispersions'].append(alpha)
 adjusted=bh([f['pValue'] for f in tested]);metrics['BH']=float(max(abs(q-f['adjustedPValue']) for q,f in zip(adjusted,tested)))
 if metrics['BH']>2e-7:failures.append(dict(failed=['BH']))
 oldcalls={f['featureIndex'] for f in old['features'] if f.get('adjustedPValue',1)<.05};newcalls={f['featureIndex'] for f in tested if f['adjustedPValue']<.05}
 result=dict(case=rec['id'],status='passed' if not failures and not unavailable else 'completed-with-failures',sourceSHA256=rec['sourceSHA256'],outputSHA256=run['outputSHA256'],checkerSHA256=sha(Path(__file__).read_bytes()),testedGenes=len(tested),oldTestedGenes=old['testedFeatures'],originalBHCalls=len(oldcalls),likelihoodRatioBHCalls=len(newcalls),gainedCallIndices=sorted(newcalls-oldcalls),lostCallIndices=sorted(oldcalls-newcalls),metrics=metrics,failures=failures,unavailable=unavailable,precisionRefinements=refinements,numpyLongDoubleBits=int(np.finfo(np.longdouble).bits))
 (d/'check.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');gzwrite(d/'gene-checks.json.gz',details)
 if rec['family']=='treatment':gzwrite(d/'r-input.json.gz',dict(groups=list(r_groups.values())))
 print(json.dumps({k:v for k,v in result.items() if k not in ['gainedCallIndices','lostCallIndices','failures','unavailable','precisionRefinements']}),flush=True)
