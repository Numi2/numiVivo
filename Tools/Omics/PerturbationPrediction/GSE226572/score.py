#!/usr/bin/env python3
"""Independent frozen-model reconstruction and all-donor/all-time outcome scoring."""
import argparse,json,importlib.metadata
from pathlib import Path
import numpy as np
from scipy.stats import t
from sklearn.kernel_ridge import KernelRidge
from download import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();src=a.study/'prediction-inputs';native=a.study/'native-prediction';a.out.mkdir(exist_ok=False)
 inf=json.loads((src/'input-freeze.json').read_text());freeze=json.loads((native/'prediction-freeze.json').read_text());assert freeze['status']=='completed-before-scoring' and not freeze['scoringStarted']
 for name,h in inf['files'].items():assert sha(src/name)==h,name
 for name,h in freeze['files'].items():assert sha(native/name)==h,name
 model=json.loads((native/'model/model.json').read_text());report=json.loads((native/'prediction/report.json').read_text());panel=json.loads((src/'panel.json').read_text());assert model['featureIDs']==report['featureIDs']==panel
 cohorts={s:json.loads((src/(s+'-cohort.json')).read_text()) for s in ['kang','query']};counts={};logs={};order={}
 for cohort in cohorts:
  with np.load(src/(cohort+'-source-counts.npz'),allow_pickle=False) as z:counts[cohort]=z['counts'];ids=z['featureIDs'].tolist()
  logs[cohort]=np.log1p(counts[cohort].astype(np.float64)/counts[cohort].sum(axis=1,keepdims=True,dtype=np.uint64)*1e6);lookup={v:i for i,v in enumerate(ids)};order[cohort]=[lookup[g] for g in panel]
 samples=cohorts['kang']['samples'];donors=sorted({s['donorID'] for s in samples});assert len(donors)==8 and donors==model['trainingDonors'];pairs=[]
 for donor in donors:
  c=[i for i,s in enumerate(samples) if s['donorID']==donor and s['condition']=='control'];r=[i for i,s in enumerate(samples) if s['donorID']==donor and s['condition']=='IFNB'];assert len(c)==len(r)==1;pairs.append((c[0],r[0]))
 c,r=np.array(pairs).T;controls=logs['kang'][c][:,order['kang']];response=logs['kang'][r][:,order['kang']]-controls;pc=counts['kang'][c][:,order['kang']]+counts['kang'][r][:,order['kang']]
 selected=np.flatnonzero((pc.sum(axis=0,dtype=np.uint64)>=10)&((pc>0).sum(axis=0)>=2));assert selected.tolist()==model['selectedFeatureIndices'];center=controls[:,selected].mean(axis=0);scale=controls[:,selected].std(axis=0);constant=np.all(controls[:,selected]==controls[0,selected],axis=0);center[constant]=controls[0,selected][constant];scale[constant|(scale==0)]=1;context=(controls[:,selected]-center)/scale/np.sqrt(len(selected));mean=response.mean(axis=0);median=np.median(response,axis=0);kernel=context@context.T;dual=np.linalg.solve(kernel+np.eye(8),response-mean);estimator=KernelRidge(alpha=1,kernel='precomputed').fit(kernel,response-mean);np.testing.assert_allclose(estimator.dual_coef_,dual,atol=1e-9,rtol=1e-8)
 model_errors={}
 for key,ref in [('contextCenters',center),('contextScales',scale),('contexts',context),('meanResponse',mean),('medianResponse',median),('dualCoefficients',dual)]:
  actual=np.asarray(model[key]);np.testing.assert_allclose(actual,ref,rtol=1e-8,atol=1e-9);model_errors[key]=float(np.abs(actual-ref).max())
 variance=response.var(axis=0,ddof=1);available=(~np.all(response==response[0],axis=0))&(variance>0);native_var=np.array([np.nan if v is None else v for v in model['donorResponseVariances']]);np.testing.assert_array_equal(np.isfinite(native_var),available);np.testing.assert_allclose(native_var[available],variance[available],atol=1e-11,rtol=1e-9);critical=float(t.ppf(.975,7));half=critical*np.sqrt(variance*(1+1/8))
 scores=[];intervals=[];comparisons=[];arrays={};qs=cohorts['query']['samples'];seen=set();assert len(report['predictions'])==3
 for prediction in report['predictions']:
  donor=prediction['group']['donorID'];assert donor not in seen and donor not in donors;seen.add(donor);q=[i for i,s in enumerate(qs) if s['donorID']==donor and s['condition']=='control'];assert len(q)==1;control=logs['query'][q[0],order['query']];np.testing.assert_allclose(prediction['control'],control,atol=1e-10,rtol=1e-10);assert prediction['libraryCounts']==int(counts['query'][q[0]].sum(dtype=np.uint64))
  query_context=(control[selected]-center)/scale/np.sqrt(len(selected));ridge=estimator.predict((query_context@context.T)[None,:])[0]+mean;refs=dict(noChange=np.zeros(len(panel)),meanResponse=mean,medianResponse=median,contextRidge=ridge);assert [e['baseline'] for e in prediction['estimates']]==list(refs)
  native_control=np.array(prediction['control']);iv=prediction['meanResponsePredictiveInterval'];assert iv['nominalCoverage']==.95 and iv['degreesOfFreedom']==7 and iv['trainingDonors']==8 and iv['unavailableFeatureIndices']==np.flatnonzero(~available).tolist();np.testing.assert_allclose(iv['studentCriticalValue'],critical,rtol=1e-11,atol=1e-11)
  bounds={};bound_refs=dict(unclippedResponseLower=mean-half,unclippedResponseUpper=mean+half,predictedTreatedLower=np.maximum(0,control+mean-half),predictedTreatedUpper=np.maximum(0,control+mean+half));maximum_bound=0.
  for key,ref in bound_refs.items():
   val=np.array([np.nan if v is None else v for v in iv[key]]);np.testing.assert_array_equal(np.isfinite(val),available);np.testing.assert_allclose(val[available],ref[available],atol=1e-9,rtol=1e-9);maximum_bound=max(maximum_bound,float(np.abs(val[available]-ref[available]).max()));bounds[key]=val;arrays[donor+'|'+key]=val
  for e in prediction['estimates']:
   method=e['baseline'];expected=np.maximum(0,control+refs[method]);actual=np.array(e['predictedTreated'])
   for key,ref in [('unclippedResponse',refs[method]),('predictedTreated',expected),('predictedResponse',expected-control)]:np.testing.assert_allclose(e[key],ref,atol=1e-8,rtol=1e-8)
   np.testing.assert_allclose(e['impliedCPMSum'],np.expm1(expected).sum(),rtol=1e-9,atol=1e-5);comparisons.append(dict(donor=donor,method=method,maximumPredictionError=float(np.abs(actual-expected).max()),maximumIntervalBoundError=maximum_bound));arrays[donor+'|'+method]=actual
  rows=[i for i,s in enumerate(qs) if s['donorID']==donor and s['condition']!='control'];assert len(rows)==6
  for i in rows:
   condition=qs[i]['condition'];hours=int(condition.removeprefix('IFNB:').removesuffix('h'));truth=logs['query'][i,order['query']];arrays[donor+'|truth|'+condition]=truth;observed=truth-native_control
   for e in prediction['estimates']:
    actual=np.array(e['predictedTreated']);change=actual-native_control;err=actual-truth;correlation=None if np.std(change)==0 or np.std(observed)==0 else float(np.corrcoef(change,observed)[0,1]);scores.append(dict(donor=donor,hours=hours,method=e['baseline'],features=len(panel),RMSE=float(np.sqrt(np.mean(err**2))),MAE=float(np.mean(np.abs(err))),responseCorrelation=correlation,clippedFeatures=int(np.count_nonzero(control+np.array(e['unclippedResponse'])<0)),panelImpliedCPM=e['impliedCPMSum'],sourceSampleIDs=cohorts['query']['sourceGroups'][i]['sampleIDs']))
   for space,y,lo,hi in [('response',observed,bounds['unclippedResponseLower'],bounds['unclippedResponseUpper']),('treated',truth,bounds['predictedTreatedLower'],bounds['predictedTreatedUpper'])]:
    v=y[available];l=lo[available];h=hi[available];intervals.append(dict(donor=donor,hours=hours,space=space,nominalCoverage=.95,availableFeatures=int(available.sum()),totalFeatures=len(panel),coverage=float(np.mean((v>=l)&(v<=h))),below=int((v<l).sum()),above=int((v>h).sum()),meanWidth=float(np.mean(h-l))))
 assert len(scores)==72 and len(intervals)==36 and len(seen)==3
 per_donor={d:{m:float(np.mean([s['RMSE'] for s in scores if s['donor']==d and s['method']==m])) for m in refs} for d in sorted(seen)};overall={m:float(np.mean([v[m] for v in per_donor.values()])) for m in refs};per_time=[dict(hours=h,donors=len({s['donor'] for s in scores if s['hours']==h}),RMSE={m:float(np.mean([s['RMSE'] for s in scores if s['hours']==h and s['method']==m])) for m in refs}) for h in sorted({s['hours'] for s in scores})]
 worse={}
 for m in refs:
  worse[m]=sum(s['RMSE']>next(v['RMSE'] for v in scores if v['donor']==s['donor'] and v['hours']==s['hours'] and v['method']=='noChange') for s in scores if s['method']==m)
 summary=dict(status='completed-with-independent-numerical-checks',trainingDonors=8,queryDonors=3,treatedDonorTimes=18,sharedFeatures=len(panel),contextFeatures=len(selected),equalDonorEqualWithinDonorTimeRMSE=overall,perDonorRMSE=per_donor,perTimeRMSE=per_time,meanGainOverNoChangePercent=100*(1-overall['meanResponse']/overall['noChange']),primaryMeanResponsePass=overall['meanResponse']<=.95*overall['noChange'],secondaryRidgePass=overall['contextRidge']<min(overall['noChange'],overall['meanResponse']),worseThanNoChangeDonorTimes=worse,intervalCoverageRange={space:[min(v['coverage'] for v in intervals if v['space']==space),max(v['coverage'] for v in intervals if v['space']==space)] for space in ['response','treated']},maximumPredictionError=max(c['maximumPredictionError'] for c in comparisons),maximumIntervalBoundError=max(c['maximumIntervalBoundError'] for c in comparisons),modelErrors=model_errors,inputFreezeSHA256=sha(src/'input-freeze.json'),predictionFreezeSHA256=sha(native/'prediction-freeze.json'),scorerSHA256=sha(Path(__file__)),packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','scikit-learn']},scope='Fixed six-hour source response across independent donor/time contexts; whole-population initial-QC RNA, no learned timing, authoritative identities, general calibration or phenotype qualification.')
 write(a.out/'summary.json',summary);write(a.out/'scores.json',scores);write(a.out/'intervals.json',intervals);write(a.out/'comparisons.json',comparisons);np.savez_compressed(a.out/'all-predictions-and-outcomes.npz',featureIDs=np.array(panel),**arrays);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
