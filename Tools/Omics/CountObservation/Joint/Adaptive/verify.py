"""Independent individual-cell likelihoods, convex dual, and query predictions.

The dual optimizes one positive variable per donor rather than native mixture
weights. Learned weights need not be unique; query arithmetic is checked for
the actual published finite-support model, not claimed uniquely identified.
"""
from pathlib import Path
import collections,gzip,hashlib,json,sys,time
import numpy as np
from scipy import sparse,optimize,special

root=Path(sys.argv[1]);origin=sys.argv[2];start=time.time()
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
inp=read(root/(origin+'-input.json.gz'));state=read(root/(origin+'-prepare-state.json'));identities=read(root/(origin+'-reference-identities.json.gz'))
assert sha(root/(origin+'-individual-cells.npz'))==inp['sourceBindings']['individualCellsSHA256']
raw=gzip.decompress((root/(origin+'-input.json.gz')).read_bytes());assert hashlib.sha256(raw).hexdigest()==state['inputSHA256']
npz=np.load(root/(origin+'-individual-cells.npz'));x=sparse.csr_matrix((npz['data'],npz['indices'],npz['indptr']),shape=tuple(npz['shape']));depth=npz['libraryCounts'];group_ids=npz['groupIndices']
genes={g['featureID']:g for g in inp['genes']};queries={q['id']:q for q in inp['queries']};raw_queries={q['id']:q for q in identities['queries']};feature_ids=npz['featureIDs'].tolist()
errors=[];max_error=0.;comparisons=0;models={};query_results={};fit_records=[];query_records=[];histograms=0
def close(actual,expected,label,tol=2e-7):
 global max_error,comparisons
 a=np.asarray(actual,float);b=np.asarray(expected,float);assert a.shape==b.shape,label
 error=float(np.max(np.abs(a-b)/np.maximum(1,np.abs(b)))) if a.size else 0.;max_error=max(max_error,error);comparisons+=int(a.size)
 if not np.isfinite(error) or error>tol:errors.append({'where':label,'maximumScaledError':error})
def histogram(y,library):
 d,inverse,n=np.unique(library,return_inverse=True,return_counts=True);total=np.bincount(inverse,weights=y).astype(np.uint64)
 return {'libraryCounts':d.tolist(),'cellsPerLibrary':n.tolist(),'geneCountsPerLibrary':total.tolist()}
for q in inp['queries']:
 raw_q=raw_queries[q['id']];assert q['control']==histogram(np.array(raw_q['counts']),np.array(raw_q['libraryCounts']));histograms+=1
 assert q['plannedLibraryCounts']==raw_q['plannedLibraryCounts'] and q['featureID']==raw_q['featureID'] and q['donorID']==raw_q['donorID']
cell_likelihoods={}
def likelihood(y,library,phi):
 e=np.asarray(library,float)/1e6;y=np.asarray(y,float);total=float(y.sum());exposure=float(e.sum())
 if total==0:mode=0.
 elif phi==0:mode=total/exposure
 else:
  def score(log_rate):
   q=special.expit(np.log(phi*e)+log_rate)
   return float(np.sum(y*(1-q)-q/phi))
  center=np.log(total/exposure);mode=float(np.exp(optimize.brentq(score,center-50,center+50,xtol=1e-13)))
 cache={}
 def ell(rate):
  rate=float(rate)
  if rate in cache:return cache[rate]
  if rate==0:value=0. if total==0 else -np.inf
  elif total==0:value=-exposure*rate if phi==0 else float(-np.log1p(phi*e*rate).sum()/phi)
  elif phi==0:value=total*np.log(rate/mode)-exposure*(rate-mode)
  else:value=total*np.log(rate/mode)-float(np.dot(y+1/phi,np.log1p(phi*e*(rate-mode)/(1+phi*e*mode))))
  cache[rate]=value;return value
 def slope(rate,scale):
  if rate==0:return -scale*exposure if total==0 else np.inf
  if phi==0:score=total-exposure*rate
  else:
   a=phi*e*rate;score=float(np.sum(y/(1+a)-(a/(1+a))/phi))
  return score*(1+scale/rate)
 return mode,ell,slope
for gene in inp['genes']:
 j=feature_ids.index(gene['featureID']);column=x[:,j].toarray().ravel()
 cs=[];ts=[]
 for pair in sorted(gene['pairs'],key=lambda p:p['donorID']):
  for condition,key,phi,out in [('control','control',gene['controlCellDispersion'],cs),('IFNB','treated',gene['treatedCellDispersion'],ts)]:
   gi=next(i for i,g in enumerate(identities['groups']) if g['donorID']==pair['donorID'] and g['conditionID']==condition);mask=group_ids==gi;y=column[mask];library=depth[mask]
   assert histogram(y,library)==pair[key];histograms+=1
   if phi is not None:out.append(likelihood(y,library,phi))
 cell_likelihoods[gene['featureID']]=(cs,ts)
def moments(cr,tr,w):
 c=np.repeat(cr,len(tr));t=np.tile(tr,len(cr));mc=float(w@c);mt=float(w@t);lc=float(w@np.log1p(c));lt=float(w@np.log1p(t))
 return {'controlMeanCPM':mc,'treatedMeanCPM':mt,'controlVarianceCPM2':float(w@((c-mc)**2)),'treatedVarianceCPM2':float(w@((t-mt)**2)),
  'covarianceCPM2':float(w@((c-mc)*(t-mt))),'responseMeanCPM':mt-mc,'responseVarianceCPM2':float(w@((t-c-(mt-mc))**2)),
  'controlMeanLog1pCPM':lc,'treatedMeanLog1pCPM':lt,'log1pResponseMean':lt-lc,'log1pResponseVariance':float(w@((np.log1p(t)-np.log1p(c)-(lt-lc))**2))}
def sampling(m,phi,libraries):
 e=np.array(libraries,float)/1e6;mean=float(e.sum()*m['treatedMeanCPM']);over=float(phi*(e@e)*(m['treatedVarianceCPM2']+m['treatedMeanCPM']**2));latent=float(e.sum()**2*m['treatedVarianceCPM2'])
 return dict(plannedCells=len(libraries),plannedLibraryCounts=sum(libraries),meanGeneCounts=mean,conditionalPoissonVariance=mean,conditionalCellOverdispersionVariance=over,latentRateVariance=latent,totalGeneCountVariance=mean+over+latent)
from fractions import Fraction
parent=Path(sys.argv[3]) if len(sys.argv)>3 else root/'parent'
old={}
with gzip.open(parent/(origin+'-native.jsonl.gz'),'rt') as f:
 for line in f:
  row=json.loads(line)
  if row['kind']=='fit' and row['modelID'].endswith('grid=65') and 'model' in row:old[row['model']['featureID']]=row['model']
certificate_records=[];available_models={};outcome=[]
with gzip.open(root/(origin+'-adaptive.jsonl.gz'),'rt') as f:
 for line in f:
  row=json.loads(line);mid=row['modelID']
  if row['kind']=='fit':
   feature=mid.split(':grid=')[0][len(origin)+1:];gene=genes[feature]
   if gene['controlCellDispersion'] is None or gene['treatedCellDispersion'] is None:
    assert row['status']=='unavailableCellDispersion' and 'model' not in row;models[mid]=None;fit_records.append({'modelID':mid,'status':row['status']});continue
   adaptive=row['model'];model=adaptive['model'];models[mid]=adaptive;cs,ts=cell_likelihoods[feature];n=len(cs)
   assert row['status']==adaptive['status'];assert bytes(model['trainingSource']['bytes']).hex()==state['inputSHA256']
   assert model['featureID']==feature and model['controlCellDispersion']==gene['controlCellDispersion'] and model['treatedCellDispersion']==gene['treatedCellDispersion']
   assert model['trainingDonorIDs']==sorted(p['donorID'] for p in gene['pairs'])
   cr=np.array(model['controlRatesCPM']);tr=np.array(model['treatedRatesCPM']);w=np.array(model['probabilities'])
   assert np.all(w>=0) and abs(w.sum()-1)<1e-10
   a=np.exp(np.array([[l[1](r) for r in cr] for l in cs])[:,:,None]+np.array([[l[1](r) for r in tr] for l in ts])[:,None,:]).reshape(n,-1)
   fitted=a@w;ll=float(np.log(fitted).sum());gap=max(0,float(np.max(a.T@(1/fitted)))/n-1)
   close(model['fittedDonorRelativeLikelihoods'],fitted,mid+'/likelihoods');close(model['relativeLogLikelihood'],ll,mid+'/objective');close(model['meanLogLikelihoodGap'],gap,mid+'/finiteGap')
   for field,value in moments(cr,tr,w).items():close(model['moments'][field],value,mid+'/'+field)
   step_ll=[r['relativeLogLikelihood'] for r in adaptive['steps']];assert all(b>=a-1e-9 for a,b in zip(step_ll,step_ll[1:]));close(step_ll[-1] if step_ll and model['status']=='convergedFiniteGridLikelihood' else ll,ll,mid+'/lastStep')
   cert=adaptive.get('certificate');upper=None;leaf_count=0
   if cert is not None:
    cscale=1/(1+gene['controlCellDispersion']*max(max(p['control']['libraryCounts']) for p in gene['pairs'])/1e6)
    tscale=1/(1+gene['treatedCellDispersion']*max(max(p['treated']['libraryCounts']) for p in gene['pairs'])/1e6)
    close(cert['controlScaleCPM'],cscale,mid+'/controlScale');close(cert['treatedScaleCPM'],tscale,mid+'/treatedScale')
    # Geometry uses declared scales and native MLE endpoint values; likelihood
    # and derivatives below use literal individual-cell counts independently.
    cscale=cert['controlScaleCPM'];tscale=cert['treatedScaleCPM']
    cm=np.array(model['controlDonorMLERatesCPM']);tm=np.array(model['treatedDonorMLERatesCPM'])
    close(cm,[v[0] for v in cs],mid+'/controlMLEs');close(tm,[v[0] for v in ts],mid+'/treatedMLEs')
    bounds=[float(np.log1p(cm.min()/cscale)),float(np.log1p(cm.max()/cscale)),float(np.log1p(tm.min()/tscale)),float(np.log1p(tm.max()/tscale))]
    caches=[{},{}]
    def axis(x,side):
     cache=caches[side]
     if x not in cache:
      likes=cs if side==0 else ts;scale=cscale if side==0 else tscale;r=scale*np.expm1(x)
      cache[x]=(r,np.array([v[1](r) for v in likes]),np.array([v[2](r,scale) for v in likes]))
     return cache[x]
    inv=1/(fitted*n);max_upper=0.;max_sample=0.;paths=sorted(b['path'] for b in cert['boxes']);assert len(set(paths))==len(paths)
    assert sum((Fraction(1,2**len(p)) for p in paths),Fraction(0))==1
    assert all(not b.startswith(a) for a,b in zip(paths,paths[1:]))
    for box in cert['boxes']:
     b=bounds.copy()
     for bit in box['path']:
      assert bit in '01';cl,ch,tl,th=b
      split=(bounds[1]>bounds[0]) and ((bounds[3]==bounds[2]) or (ch-cl)/(bounds[1]-bounds[0]) >= (th-tl)/(bounds[3]-bounds[2]))
      i=0 if split else 2;midpoint=(b[i]+b[i+1])/2;b[i if bit=='1' else i+1]=midpoint
     actual=[box[k] for k in ['controlLow','controlHigh','treatedLow','treatedHigh']];close(actual,b,mid+'/partition/'+box['path'],tol=1e-12)
     x0,x1,y0,y1=actual;x=(x0+x1)/2;y=(y0+y1)/2
     ar0,a0,_=axis(x0,0);ar1,a1,_=axis(x1,0);br0,b0,_=axis(y0,1);br1,b1,_=axis(y1,1);_,am,sa=axis(x,0);_,bm,sb=axis(y,1)
     ind=float(np.sum(np.exp(np.where((cm>=ar0)&(cm<=ar1),0,np.maximum(a0,a1))+np.where((tm>=br0)&(tm<=br1),0,np.maximum(b0,b1)))*inv))
     tangent=0.
     for xx in [x0,x1]:
      for yy in [y0,y1]:
       dx=0 if xx==x else sa*(xx-x);dy=0 if yy==y else sb*(yy-y)
       with np.errstate(over='ignore',invalid='ignore'):value=float(np.sum(np.exp(am+bm+dx+dy)*inv))
       tangent=max(tangent,value) if not np.isnan(value) else np.inf
     raw=min(ind,tangent);bound=raw+1e-9*(1+raw);close(box['meanDirectionalUpperBound'],bound,mid+'/boxBound/'+box['path']);max_upper=max(max_upper,bound)
     sample=float(np.sum(np.exp(am+bm)*inv));assert sample<=bound+1e-8*max(1,bound);max_sample=max(max_sample,sample)
    witness=float(sum(np.exp(cs[d][1](cert['bestControlRateCPM'])+ts[d][1](cert['bestTreatedRateCPM']))*inv[d] for d in range(n)))
    close(cert['maximumMeanDirectionalLowerBound'],witness,mid+'/witness');close(cert['maximumMeanDirectionalUpperBound'],max_upper,mid+'/globalBound')
    upper=max_upper;leaf_count=len(paths)
    if adaptive['status']=='boundedContinuousLikelihood':assert upper<=1+adaptive['plan']['meanLogLikelihoodGapTolerance']+2e-8
    certificate_records.append({'modelID':mid,'leaves':leaf_count,'completePrefixFreePartition':True,'independentUpper':upper,'independentWitness':witness,'maximumLeafCenterScore':max_sample})
   previous=old[feature];ll_gain=(ll-previous['relativeLogLikelihood'])/n
   if adaptive['status']=='boundedContinuousLikelihood':assert ll_gain>=-adaptive['plan']['meanLogLikelihoodGapTolerance']-2e-7
   fit_records.append({'modelID':mid,'status':row['status'],'leaves':leaf_count,'independentUpper':upper,'meanLogLikelihoodGainOverGrid65':ll_gain,'supportAdditions':len(adaptive['steps'])-1,'supportShape':[len(cr),len(tr)],'finiteGridGap':gap})
  else:
   q=queries[row['queryID']];adaptive=models[mid]
   if adaptive is None or adaptive['status']!='boundedContinuousLikelihood':
    expected='unavailableCellDispersion' if adaptive is None else adaptive['status'];assert row['status']==expected and 'prediction' not in row;query_records.append({'modelID':mid,'queryID':q['id'],'status':row['status']});continue
   assert row['status']=='conditionalPrediction';model=adaptive['model'];pr=row['prediction'];rawq=raw_queries[q['id']];_,ell,_=likelihood(np.array(rawq['counts']),np.array(rawq['libraryCounts']),model['controlCellDispersion'])
   cr=np.array(model['controlRatesCPM']);tr=np.array(model['treatedRatesCPM']);w=np.array(model['probabilities']);logw=np.full(w.shape,-np.inf);positive=w>0;logw[positive]=np.log(w[positive])+np.repeat([ell(r) for r in cr],len(tr))[positive];post=np.exp(logw-special.logsumexp(logw));close(pr['probabilities'],post,mid+'/'+q['id']+'/queryWeights')
   assert bytes(pr['querySource']['bytes']).hex()==state['inputSHA256'] and pr['featureID']==q['featureID'] and pr['queryDonorID']==q['donorID'] and pr['treatedConditionID']=='IFNB'
   supported=np.isfinite(logw);coordinates=np.tile(tr,len(cr));degenerate=len(np.unique(coordinates[supported]))==1
   assert pr['degenerateTreatedRateDistribution']==degenerate
   assert pr['underflowedPosteriorComponents']==int((supported&(np.array(pr['probabilities'])==0)).sum())
   assert pr['numericallyDegenerateTreatedRateMoments']==(not degenerate and pr['moments']['treatedVarianceCPM2']==0)
   if degenerate:assert pr['moments']['treatedVarianceCPM2']==0
   mm=moments(cr,tr,post);sample=sampling(mm,model['treatedCellDispersion'],q['plannedLibraryCounts'])
   for k,v in mm.items():close(pr['moments'][k],v,mid+'/'+q['id']+'/'+k)
   for k,v in sample.items():close(pr['plannedTreatedCountMoments'][k],v,mid+'/'+q['id']+'/'+k)
   query_records.append({'modelID':mid,'queryID':q['id'],'status':row['status'],'moments':mm,'sample':sample})
assert len(fit_records)==16 and len(query_records)==992
result={'status':'passed-independent-numerics' if not errors else 'failed-independent-numerics','origin':origin,'histogramsChecked':histograms,'numericalValuesCompared':comparisons,'maximumScaledError':max_error,'errors':errors,'fits':fit_records,'certificates':certificate_records,'fitStates':dict(collections.Counter(r['status'] for r in fit_records)),'queryStates':dict(collections.Counter(r['status'] for r in query_records)),'seconds':time.time()-start,'qualification':'Analytic support bounds reconstructed in FP64 from individual cells, not directed-rounding interval certification or biological validation.'}
(root/(origin+'-verification.json')).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
with (root/(origin+'-query-reference.json.gz')).open('wb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(json.dumps(query_records,sort_keys=True,separators=(',',':'),allow_nan=False).encode())
print(json.dumps({k:v for k,v in result.items() if k not in ['fits','certificates']}));sys.exit(bool(errors))
