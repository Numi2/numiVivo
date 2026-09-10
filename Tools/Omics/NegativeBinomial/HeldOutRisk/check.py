#!/usr/bin/env python3
"""Independent training isolation, conditional fit and held-out score checks."""
import argparse,gzip,hashlib,json,shlex,subprocess
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import nbinom,norm

def load(p):return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dense(b):
 m=b['matrix'];return csr_matrix((np.array(m['counts'],dtype=np.int64),m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount'])).toarray()
TOL=dict(scaledScore=1e-7,standardDeviation=2e-6,aggregateLogScore=2e-6,meanSquaredLog1pError=2e-10)

def check(root,rec,original,source_counts,model_root=None,model_host='macmini'):
 d=root/rec['id'];training=load(d/'training.json.gz');test=load(d/'test.json.gz');scores=load(d/'scores.json.gz')
 if model_root:
  assert len(rec['id'])==3 and rec['id'].isdigit()
  model_bytes=subprocess.check_output(['ssh',model_host,'cat '+shlex.quote(str(Path(model_root)/rec['id']/'model.json.gz'))])
 else:model_bytes=(d/'model.json.gz').read_bytes()
 model_sha=hashlib.sha256(model_bytes).hexdigest();model=json.loads(gzip.decompress(model_bytes))
 run=next(r for r in load(d/'runs.json') if r['stage']=='fit')
 assert model_sha==run['outputSHA256'] and hashlib.sha256(gzip.decompress(model_bytes)).hexdigest()==run['logicalSHA256']
 assert sha(d/'training.json.gz')==rec['trainingSHA256'] and sha(d/'test.json.gz')==rec['testSHA256']
 counts=dense(training['bulk']);outcomes=np.array(test['counts'],dtype=np.int64)
 np.testing.assert_array_equal(counts,source_counts[rec['trainingSourceRows']]);np.testing.assert_array_equal(outcomes,source_counts[rec['testSourceRows']])
 assert training['bulk']['groups']==[original['pseudobulk']['groups'][i] for i in rec['trainingSourceRows']]
 expected_request=dict(training['request'])
 expected_request['negativeBinomialOptions']={**dict(minimumTrendGenes=20,minimumPriorVariance=.25,outlierStandardDeviations=2),**expected_request['negativeBinomialOptions']}
 assert model['request']==expected_request and model['request']['includedDonorIDs']==rec['trainingAnimals']
 assert set(rec['trainingAnimals']).isdisjoint(rec['testAnimals'])
 design=model['design'];rows=design['sourcePseudobulkIndices'];counts=counts[rows]
 assert design['controlReplicates']==design['treatmentReplicates']==3
 assert design['residualDegreesOfFreedom']==4
 x=np.array(design['rows']);np.testing.assert_array_equal(x[:,0],1)
 np.testing.assert_array_equal(x[:,1],[g['condition']=='LPS' for g in design['observations']])
 assert x.shape==(6,2) and design['contrast']==[0,1]
 libraries=counts.sum(axis=1);center=np.log(libraries).mean();factors=np.exp(np.log(libraries)-center)
 np.testing.assert_allclose(factors,design['sizeFactorValues'],rtol=2e-14,atol=1e-15)
 assert abs(center-model['trainingLogLibraryCenter'])<2e-14
 genes=model['genes'];assert [g['featureIndex'] for g in genes]==list(range(counts.shape[1]))
 filtered=(counts.sum(axis=0)<10)|((counts>0).sum(axis=0)<3)
 np.testing.assert_array_equal(filtered,[g['status']=='filteredLowExpression' for g in genes])
 mean=(counts/factors[:,None]).mean(axis=0);np.testing.assert_allclose(mean,[g['mean'] for g in genes],rtol=2e-14,atol=2e-12)
 tested=[g for g in genes if g['status']=='tested'];idx=[g['featureIndex'] for g in tested]
 y=counts[:,idx].T;alpha=np.array([g['dispersion'] for g in tested]);ref_genes=[g for g in tested if abs(g['mle']['coefficients'][1]/np.log(2))<10]
 prior=model.get('prior');assert prior and not model.get('priorError')
 assert prior['eligibleFeatureIndices']==idx and prior['referenceFeatureIndices']==[g['featureIndex'] for g in ref_genes]
 weights=np.array([1/(1/g['mean']+g['trendDispersion']) for g in ref_genes]);effects=np.abs([g['mle']['coefficients'][1]/np.log(2) for g in ref_genes])
 order=np.argsort(effects,kind='stable');cum=np.cumsum(weights[order]*len(weights)/weights.sum());rank=1+(cum[-1]-1)*.95;low=np.floor(rank);high=min(low+1,cum[-1]);fraction=rank-low
 quantile=(1-fraction)*effects[order[min(np.searchsorted(cum,low),len(order)-1)]]+fraction*effects[order[min(np.searchsorted(cum,high),len(order)-1)]]
 sd=max(.01,quantile/norm.ppf(.975));assert abs(sd-prior['priorStandardDeviationLog2'])<2e-12
 errors=[];metrics=[];arm=x[:,1]
 for method in ['mle','fixed','empirical']:
  if any(g.get(method) is None or not g[method]['converged'] for g in tested):continue
  beta=np.array([g[method]['coefficients'] for g in tested]);mu=np.exp(np.log(factors)[None,:]+beta@x.T)
  precision=0 if method=='mle' else 1/(np.log(2)*(1 if method=='fixed' else sd))**2
  raw=(y-mu)/(1+alpha[:,None]*mu);score=raw@x;score[:,1]-=precision*beta[:,1]
  w=mu/(1+alpha[:,None]*mu) if method=='mle' else (1+alpha[:,None]*y)*mu/(1+alpha[:,None]*mu)**2
  h00=w.sum(axis=1);h01=w@arm;h11=w@(arm*arm)+precision
  scaled=np.max(np.abs(score)/np.sqrt(np.stack([h00,h11],axis=1)),axis=1)
  ref_sd=np.sqrt(h00/(h00*h11-h01*h01));actual_sd=np.array([g[method]['standardDeviation'] for g in tested])
  values=dict(scaledScore=scaled,standardDeviation=np.abs(ref_sd-actual_sd))
  for key,values_array in values.items():
   worst=int(np.argmax(values_array));metrics.append(dict(method=method,metric=key,maximum=float(values_array[worst]),featureIndex=idx[worst]))
   for i in np.flatnonzero(~np.isfinite(values_array)|(values_array>TOL[key])):errors.append(dict(method=method,metric=key,featureIndex=idx[i],error=float(values_array[i])))
 assert scores['testedTrainingGenes']==len(tested)
 expected_unavailable=[g['featureIndex'] for g in tested if any(not g.get(k,{}).get('converged',False) for k in ['mle','fixed','empirical'])]
 assert scores['unavailableMAPFeatureIndices']==expected_unavailable
 assert scores['primaryAvailable']==(bool(tested) and not expected_unavailable)
 ref_scores=[]
 if scores['primaryAvailable']:
  for i,sample in enumerate(test['sampleIDs']):
   library=int(outcomes[i].sum());factor=np.exp(np.log(library)-center);arm_value=int(test['conditions'][i]=='LPS')
   for method,native_name in [('mle','MLE'),('fixed','fixed'),('empirical','empirical')]:
    beta=np.array([g[method]['coefficients'] for g in tested]);prediction=np.exp(beta[:,0]+arm_value*beta[:,1]);mu=factor*prediction
    masses=nbinom.logpmf(outcomes[i,idx],1/alpha,1/(1+alpha*mu));value=float(masses.mean())
    squared=float(np.mean((np.log1p(outcomes[i,idx]/factor)-np.log1p(prediction))**2))
    native=next(o for o in scores['observations'] if o['sampleID']==sample and o['method']==native_name)
    assert native['genes']==len(tested) and native['libraryCounts']==library
    np.testing.assert_allclose(native['sizeFactor'],factor,rtol=2e-14)
    differences=dict(aggregateLogScore=abs(value-native['meanLogScore']),meanSquaredLog1pError=abs(squared-native['meanSquaredLog1pError']))
    for key,error in differences.items():
     metrics.append(dict(method=method,metric=key,maximum=error,sampleID=sample))
     if not np.isfinite(error) or error>TOL[key]:errors.append(dict(method=method,metric=key,sampleID=sample,error=error))
    ref_scores.append(dict(sampleID=sample,condition=test['conditions'][i],method=native_name,genes=len(tested),meanLogScore=value,meanSquaredLog1pError=squared))
 result=dict(case=rec['id'],cellGroup=rec['cellGroup'],status='passed' if not errors else 'failed',errors=errors,metrics=metrics,
  testedGenes=len(tested),conditionalFits=3*len(tested),trainingStatusCounts=dict(Counter(g['status'] for g in genes)),
  priorSDLog2=sd,primaryAvailable=scores['primaryAvailable'],unavailableMAPFeatureIndices=expected_unavailable,observations=ref_scores,
  sourceCountsAndIsolationVerified=True,modelSHA256=model_sha,scoresSHA256=sha(d/'scores.json.gz'),checkerSHA256=sha(Path(__file__)))
 (d/'check.json').write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n');return result

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--case');p.add_argument('--model-root');p.add_argument('--model-host',default='macmini');a=p.parse_args()
 protocol=load(a.root/'protocol.json');source=Path(protocol['sourceReport']);assert sha(source)==protocol['sourceReportSHA256'];original=load(source);counts=dense(original['pseudobulk'])
 results=[]
 for rec in protocol['cases']:
  if a.case is None or rec['id']==a.case:
   result=check(a.root,rec,original,counts,a.model_root,a.model_host);results.append(result);print(json.dumps(dict(case=rec['id'],status=result['status'],testedGenes=result['testedGenes'],errors=len(result['errors']))),flush=True)
 assert results and all(r['status']=='passed' for r in results)
