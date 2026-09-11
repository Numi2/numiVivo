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
 return mode,ell
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
with gzip.open(root/(origin+'-native.jsonl.gz'),'rt') as f:
 for line in f:
  row=json.loads(line);mid=row['modelID']
  if row['kind']=='fit':
   feature=mid.split(':grid=')[0][len(origin)+1:];gene=genes[feature]
   if gene['controlCellDispersion'] is None or gene['treatedCellDispersion'] is None:
    assert row['status']=='unavailableCellDispersion' and 'model' not in row;models[mid]=None;fit_records.append({'modelID':mid,'status':row['status']});continue
   model=row['model'];models[mid]=model;assert row['status']==model['status'];assert model['featureID']==feature
   assert bytes(model['trainingSource']['bytes']).hex()==state['inputSHA256']
   assert model['trainingDonorIDs']==sorted(p['donorID'] for p in gene['pairs'])
   assert model['controlCellDispersion']==gene['controlCellDispersion'] and model['treatedCellDispersion']==gene['treatedCellDispersion']
   assert model['controlCells']==[sum(p['control']['cellsPerLibrary']) for p in sorted(gene['pairs'],key=lambda p:p['donorID'])]
   assert model['treatedCells']==[sum(p['treated']['cellsPerLibrary']) for p in sorted(gene['pairs'],key=lambda p:p['donorID'])]
   cs,ts=cell_likelihoods[feature];close(model['controlDonorMLERatesCPM'],[l[0] for l in cs],mid+'/controlModes');close(model['treatedDonorMLERatesCPM'],[l[0] for l in ts],mid+'/treatedModes')
   cr=np.array(model['controlRatesCPM']);tr=np.array(model['treatedRatesCPM']);w=np.array(model['probabilities']);n=len(cs)
   for rates,mles in [(cr,model['controlDonorMLERatesCPM']),(tr,model['treatedDonorMLERatesCPM'])]:
    expected=set(mles);low=min(mles);high=max(mles);size=model['plan']['gridPointsPerAxis']
    if low!=high:
     for i in range(1,size-1):expected.add(float(np.expm1(np.log1p(low)+(np.log1p(high)-np.log1p(low))*i/(size-1))))
    close(rates,sorted(expected),mid+'/grid')
   a=np.exp(np.array([[l[1](r) for r in cr] for l in cs])[:,:,None]+np.array([[l[1](r) for r in tr] for l in ts])[:,None,:]).reshape(n,-1)
   assert np.all(w>=0) and abs(w.sum()-1)<1e-10
   fitted=a@w;assert np.all(fitted>0);ll=float(np.log(fitted).sum());gap=max(0,float(np.max(a.T@(1/fitted)))/n-1)
   close(model['fittedDonorRelativeLikelihoods'],fitted,mid+'/donorLikelihoods');close(model['relativeLogLikelihood'],ll,mid+'/logLikelihood');close(model['meanLogLikelihoodGap'],gap,mid+'/certificate')
   if model['status']=='convergedFiniteGridLikelihood':assert gap<=model['plan']['meanLogLikelihoodGapTolerance']+2e-9,(mid,gap)
   elif model['status']=='iterationLimit':assert model['iterations']==model['plan']['maximumIterations']
   # Independent dual, only n variables; no reconstruction of native weights.
   dual=optimize.minimize(lambda v:-np.log(v*n).sum(),np.ones(n)/n,jac=lambda v:-1/v,method='SLSQP',bounds=[(1/n,1)]*n,
       constraints=[{'type':'ineq','fun':lambda v:1-a.T@v,'jac':lambda v:-a.T}],options={'ftol':1e-12,'maxiter':1000})
   u=dual.x*n;upper=float(np.max(a.T@u)-np.log(u).sum()-n);dual_gap=(upper-ll)/n
   if model['status']=='convergedFiniteGridLikelihood' and (dual_gap < -2e-7 or dual_gap>2e-7):errors.append({'where':mid+'/independentDualGap','value':dual_gap})
   mm=moments(cr,tr,w)
   for key,value in mm.items():close(model['moments'][key],value,mid+'/'+key)
   fit_records.append({'modelID':mid,'status':row['status'],'gridShape':[len(cr),len(tr)],'positiveWeights':int((w>0).sum()),'iterations':model['iterations'],'nativeMeanGap':model['meanLogLikelihoodGap'],
    'independentRelativeLogLikelihood':ll,'independentDualUpperBound':upper,'independentMeanDualGap':dual_gap,'dualVariables':u.tolist(),'independentDonorRelativeLikelihoods':fitted.tolist(),'dualOptimizerSuccess':bool(dual.success),'dualOptimizerMessage':str(dual.message),'moments':mm})
  else:
   q=queries[row['queryID']];model=models[mid]
   if model is None or model['status']!='convergedFiniteGridLikelihood':
    expected='unavailableCellDispersion' if model is None else model['status'];assert row['status']==expected and 'prediction' not in row
    query_records.append({'modelID':mid,'queryID':q['id'],'status':row['status']});continue
   assert row['status']=='conditionalPrediction',(mid,q['id'],row)
   prediction=row['prediction'];raw_q=raw_queries[q['id']];_,ell=likelihood(np.array(raw_q['counts']),np.array(raw_q['libraryCounts']),model['controlCellDispersion'])
   cr=np.array(model['controlRatesCPM']);tr=np.array(model['treatedRatesCPM']);w=np.array(model['probabilities']);logw=np.full(w.shape,-np.inf);positive=w>0
   logw[positive]=np.log(w[positive])+np.repeat([ell(r) for r in cr],len(tr))[positive];post=np.exp(logw-special.logsumexp(logw))
   close(prediction['probabilities'],post,mid+'/'+q['id']+'/posteriorWeights');mm=moments(cr,tr,post);sample=sampling(mm,model['treatedCellDispersion'],q['plannedLibraryCounts'])
   for key,value in mm.items():close(prediction['moments'][key],value,mid+'/'+q['id']+'/'+key)
   for key,value in sample.items():close(prediction['plannedTreatedCountMoments'][key],value,mid+'/'+q['id']+'/'+key)
   assert bytes(prediction['querySource']['bytes']).hex()==state['inputSHA256'] and prediction['featureID']==q['featureID'] and prediction['queryDonorID']==q['donorID'] and prediction['treatedConditionID']=='IFNB'
   supported=np.isfinite(logw);coordinates=np.tile(tr,len(cr));degenerate=len(np.unique(coordinates[supported]))==1
   assert prediction['degenerateTreatedRateDistribution']==degenerate
   native_weights=np.array(prediction['probabilities']);assert prediction['underflowedPosteriorComponents']==int((supported&(native_weights==0)).sum())
   assert prediction['numericallyDegenerateTreatedRateMoments']==(not degenerate and prediction['moments']['treatedVarianceCPM2']==0)
   if degenerate:assert prediction['moments']['treatedVarianceCPM2']==0
   record={'modelID':mid,'queryID':q['id'],'status':row['status'],'moments':mm,'sample':sample,'degenerate':degenerate,'underflowedPosteriorComponents':prediction['underflowedPosteriorComponents'],'numericallyDegenerate':prediction['numericallyDegenerateTreatedRateMoments']};query_records.append(record);query_results[(model['featureID'],model['plan']['gridPointsPerAxis'],q['id'])]=record
assert len(fit_records)==16*len(inp['gridSizes']) and len(query_records)==992*len(inp['gridSizes'])
refine=read(root/'refinement-protocol.json');refinement=[]
for gene in inp['genes']:
 feature=gene['featureID'];fine=models[f'{origin}:{feature}:grid=65'];coarse=models[f'{origin}:{feature}:grid=33']
 if fine is None or coarse is None:refinement.append({'featureID':feature,'status':'unavailableCellDispersion'});continue
 if fine['status']!='convergedFiniteGridLikelihood' or coarse['status']!='convergedFiniteGridLikelihood':refinement.append({'featureID':feature,'status':'fitNotConverged'});continue
 llchange=(fine['relativeLogLikelihood']-coarse['relativeLogLikelihood'])/len(fine['trainingDonorIDs']);maximum=0.;worst=None
 for q in inp['queries']:
  if q['featureID']!=feature:continue
  f=query_results[(feature,65,q['id'])];c=query_results[(feature,33,q['id'])]
  for field in refine['fields']:
   a=f['sample'][field] if field in f['sample'] else f['moments'][field];b=c['sample'][field] if field in c['sample'] else c['moments'][field]
   error=abs(a-b)/max(1,abs(a))
   if error>maximum:maximum=error;worst={'queryID':q['id'],'field':field,'fine':a,'coarse':b}
 ok=-2e-8<=llchange<=refine['maximumMeanDonorLogLikelihoodChange'] and maximum<=refine['maximumScaledQueryMomentChange']
 refinement.append({'featureID':feature,'status':'passed' if ok else 'failedGridRefinement','meanDonorLogLikelihoodChange':llchange,'maximumScaledQueryMomentChange':maximum,'worst':worst})
result={'status':'passed-independent-numerics' if not errors else 'failed-independent-numerics','origin':origin,'histogramsChecked':histograms,'numericalValuesCompared':comparisons,'maximumScaledError':max_error,'errors':errors,
 'fits':fit_records,'fitStates':dict(collections.Counter(f['status'] for f in fit_records)),'queryStates':dict(collections.Counter(q['status'] for q in query_records)),
 'degeneratePredictions':sum(q.get('degenerate',False) for q in query_records),'refinement':refinement,'seconds':time.time()-start,
 'qualification':'Finite-support conditional likelihood and arithmetic only; parameter uncertainty, continuous mixing law, full transcriptome scaling and biological prediction calibration remain open.'}
with (root/(origin+'-query-reference.json.gz')).open('wb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(json.dumps(query_records,sort_keys=True,separators=(',',':'),allow_nan=False).encode())
(root/(origin+'-verification.json')).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['fits','refinement']}));sys.exit(bool(errors))
