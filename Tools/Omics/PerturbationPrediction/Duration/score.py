#!/usr/bin/env python3
"""Independent whole-coordinate arithmetic checks and every held-out outcome."""
import argparse,importlib.metadata,json
from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import t
from scipy import sparse
from common import sha,write,verify_files,hours

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();src=a.study/'inputs';native=a.study/'native';a.out.mkdir(exist_ok=False)
 f=json.loads((src/'input-freeze.json').read_text());frozen=json.loads((native/'prediction-freeze.json').read_text());assert frozen['status']=='completed-before-duration-scoring';verify_files(src,f['files']);verify_files(native,frozen['files'])
 cohort=json.loads((src/'query-cohort.json').read_text());samples=cohort['samples'];panel=json.loads((src/'panel.json').read_text());folds=json.loads((src/'folds.json').read_text())
 with np.load(src/'query-source-counts.npz',allow_pickle=False) as z:counts=z['counts'];features=z['featureIDs'].tolist()
 logs=np.log1p(counts.astype(np.float64)/counts.sum(axis=1,keepdims=True,dtype=np.uint64)*1e6);lookup={v:i for i,v in enumerate(features)};order=np.array([lookup[g] for g in panel]);m=len(panel)
 scores=[];intervals=[];comparisons=[];arrays={};models=[];critical=float(t.ppf(.975,1));seen=set()
 for fold in folds:
  fid=fold['id'];donor=fold['heldOutDonor'];seen.add(donor);base=native/fid;model=json.loads((base/'model/model.json').read_text());report=json.loads((base/'prediction/report.json').read_text());bulk=json.loads((base/'model/training/report.json').read_text())['pseudobulk'];matrix=bulk['matrix']
  native_counts=sparse.csr_matrix((np.array(matrix['counts'],dtype=np.uint64),matrix['featureIndices'],matrix['rowOffsets']),shape=(matrix['cellCount'],matrix['featureCount'])).toarray()
  assert bulk['featureIDs']==features and model['featureIDs']==report['featureIDs']==panel and len(bulk['groups'])==14
  for row,g in enumerate(bulk['groups']):
   assert len(g['sampleIDs'])==1;idx=next(i for i in fold['trainingRows'] if samples[i]['id']==g['sampleIDs'][0]);np.testing.assert_array_equal(native_counts[row],counts[idx]);assert g['donorID']==samples[idx]['donorID'] and g['condition']==samples[idx]['condition']
  donors=sorted({samples[i]['donorID'] for i in fold['trainingRows']});assert len(donors)==2 and donor not in donors and [c['donorID'] for c in model['curves']]==donors
  curves=[];controlrows=[];perdonorcounts=[];invariant=[];modelerrors={}
  for d,nativecurve in zip(donors,model['curves']):
   rows=[i for i in fold['trainingRows'] if samples[i]['donorID']==d];c=[i for i in rows if samples[i]['condition']=='control'];r=sorted([i for i in rows if samples[i]['condition']!='control'],key=lambda i:hours(samples[i]['condition']));assert len(c)==1 and len(r)==6
   curve=np.vstack([np.zeros(m),logs[r][:,order]-logs[c[0],order]]);times=np.array([0]+[hours(samples[i]['condition']) for i in r]);np.testing.assert_array_equal(nativecurve['hours'],times);np.testing.assert_allclose(nativecurve['responses'],curve,atol=1e-8,rtol=1e-8)
   curves.append(interp1d(np.log1p(times),curve,axis=0,kind='linear',bounds_error=True));controlrows.append(c[0]);perdonorcounts.append(counts[rows][:,order].sum(axis=0,dtype=np.uint64));invariant.append(curve[1:].mean(axis=0))
  pc=np.array(perdonorcounts);selected=np.flatnonzero((pc.sum(axis=0,dtype=np.uint64)>=10)&((pc>0).sum(axis=0)>=2));assert selected.tolist()==model['selectedFeatureIndices'];controls=logs[controlrows][:,order];center=controls[:,selected].mean(axis=0);scale=controls[:,selected].std(axis=0);constant=np.all(controls[:,selected]==controls[0,selected],axis=0);center[constant]=controls[0,selected][constant];scale[constant|(scale==0)]=1;context=(controls[:,selected]-center)/scale/np.sqrt(len(selected));system=context@context.T+np.eye(2);inverse=np.linalg.solve(system,np.eye(2))
  for key,ref in [('contextCenters',center),('contextScales',scale),('contexts',context),('ridgeInverse',inverse)]:
   actual=np.array(model[key]);np.testing.assert_allclose(actual,ref,atol=1e-8,rtol=1e-8);modelerrors[key]=float(np.max(np.abs(actual-ref)))
  assert model['maximumSupportedHours']==36;models.append(dict(fold=fid,trainingDonors=donors,selectedFeatures=len(selected),maximumInverseResidual=model['maximumSolveResidual'],errors=modelerrors,nativeSourceCountsExact=True))
  controlrow=fold['queryRows'][0];control=logs[controlrow,order];querycontext=(control[selected]-center)/scale/np.sqrt(len(selected));invariantmean=np.mean(invariant,axis=0);matched=np.maximum(0,control+invariantmean)
  assert len(report['predictions'])==6 and [v['hours'] for v in report['predictions']]==fold['hours']
  for item in report['predictions']:
   h=item['hours'];prediction=item['prediction'];assert prediction['group']['donorID']==donor and prediction['group']['condition']=='control';np.testing.assert_allclose(prediction['control'],control,atol=1e-10,rtol=1e-10);assert prediction['libraryCounts']==int(counts[controlrow].sum(dtype=np.uint64))
   responses=np.array([curve(np.log1p(h)) for curve in curves]);mean=responses.mean(axis=0);median=np.median(responses,axis=0);ridge=mean+(querycontext@context.T)@np.linalg.solve(system,responses-mean)
   references=dict(noChange=np.zeros(m),durationMean=mean,durationMedian=median,durationContextRidge=ridge);assert [e['baseline'] for e in prediction['estimates']]==list(references)
   rows=[i for i in fold['outcomeRows'] if hours(samples[i]['condition'])==h];assert len(rows)==1;truthrow=rows[0];truth=logs[truthrow,order];observed=truth-control;key=fid+'|'+str(h);arrays[key+'|truth']=truth;arrays[key+'|control']=control
   allestimates=[]
   for e in prediction['estimates']:
    method=e['baseline'];reference=references[method];expected=np.maximum(0,control+reference);actual=np.array(e['predictedTreated']);maxerror=0.
    for name,ref in [('unclippedResponse',reference),('predictedTreated',expected),('predictedResponse',expected-control)]:
     np.testing.assert_allclose(e[name],ref,atol=1e-8,rtol=1e-8);maxerror=max(maxerror,float(np.max(np.abs(np.array(e[name])-ref))))
    np.testing.assert_allclose(e['impliedCPMSum'],np.expm1(expected).sum(),atol=1e-5,rtol=1e-9);comparisons.append(dict(fold=fid,hours=h,method=method,maximumPointError=maxerror));allestimates.append((method,actual,np.array(e['unclippedResponse'])));arrays[key+'|'+method]=actual
   allestimates.append(('matchedInvariantMean',matched,invariantmean));arrays[key+'|matchedInvariantMean']=matched
   for method,actual,response in allestimates:
    error=actual-truth;applied=actual-control;correlation=None if np.std(applied)==0 or np.std(observed)==0 else float(np.corrcoef(applied,observed)[0,1]);scores.append(dict(donor=donor,hours=h,method=method,features=m,RMSE=float(np.sqrt(np.mean(error**2))),MAE=float(np.mean(np.abs(error))),responseCorrelation=correlation,clippedFeatures=int(np.count_nonzero(control+response<0)),sourceSampleIDs=cohort['sourceGroups'][truthrow]['sampleIDs']))
   variance=responses.var(axis=0,ddof=1);available=(~np.all(responses==responses[0],axis=0))&(variance>0);half=critical*np.sqrt(variance*1.5);iv=prediction['meanResponsePredictiveInterval'];assert iv['nominalCoverage']==.95 and iv['trainingDonors']==2 and iv['degreesOfFreedom']==1 and iv['unavailableFeatureIndices']==np.flatnonzero(~available).tolist();np.testing.assert_allclose(iv['studentCriticalValue'],critical,atol=1e-9,rtol=1e-9)
   bounds={};bounderrors=[]
   for name,ref in [('unclippedResponseLower',mean-half),('unclippedResponseUpper',mean+half),('predictedTreatedLower',np.maximum(0,control+mean-half)),('predictedTreatedUpper',np.maximum(0,control+mean+half))]:
    actual=np.array([np.nan if v is None else v for v in iv[name]]);np.testing.assert_array_equal(np.isfinite(actual),available);np.testing.assert_allclose(actual[available],ref[available],atol=1e-8,rtol=1e-8);bounderrors.append(float(np.max(np.abs(actual[available]-ref[available]))) if np.any(available) else 0.);bounds[name]=actual;arrays[key+'|'+name]=actual
   for space,y,lo,hi in [('response',observed,bounds['unclippedResponseLower'],bounds['unclippedResponseUpper']),('treated',truth,bounds['predictedTreatedLower'],bounds['predictedTreatedUpper'])]:
    v=y[available];l=lo[available];u=hi[available];intervals.append(dict(donor=donor,hours=h,space=space,nominalCoverage=.95,availableFeatures=int(available.sum()),unavailableFeatures=int((~available).sum()),coverage=float(np.mean((v>=l)&(v<=u))),meanWidth=float(np.mean(u-l)),medianWidth=float(np.median(u-l)),maximumBoundError=max(bounderrors)))
 assert len(seen)==3 and len(scores)==90 and len(intervals)==36 and len(comparisons)==72
 methods=list(references)+['matchedInvariantMean'];perdonor={d:{method:float(np.mean([s['RMSE'] for s in scores if s['donor']==d and s['method']==method])) for method in methods} for d in sorted(seen)};overall={method:float(np.mean([v[method] for v in perdonor.values()])) for method in methods};perhour=[dict(hours=h,donors=len({s['donor'] for s in scores if s['hours']==h}),RMSE={method:float(np.mean([s['RMSE'] for s in scores if s['hours']==h and s['method']==method])) for method in methods}) for h in sorted({s['hours'] for s in scores})]
 def worse(method,baseline):
  return sum(s['RMSE']>next(r['RMSE'] for r in scores if r['donor']==s['donor'] and r['hours']==s['hours'] and r['method']==baseline) for s in scores if s['method']==method)
 summary=dict(status='completed-with-independent-numerical-checks',evidence='retrospective-model-development-not-independent-validation',folds=3,trainingDonorsPerFold=2,heldOutDonors=3,treatedDonorTimes=18,sourceCells=cohort['sourceCells'],sharedFeatures=m,equalDonorEqualTimeRMSE=overall,perDonorRMSE=perdonor,perTimeRMSE=perhour,meanImprovementPercent={b:100*(1-overall['durationMean']/overall[b]) for b in ['noChange','matchedInvariantMean']},primaryDurationMeanPass=all(overall['durationMean']<=.95*overall[b] for b in ['noChange','matchedInvariantMean']),secondaryRidgePass=overall['durationContextRidge']<min(overall['noChange'],overall['durationMean']),worseDonorTimes={method:{b:worse(method,b) for b in ['noChange','matchedInvariantMean']} for method in methods},intervalCoverageRange={space:[min(v['coverage'] for v in intervals if v['space']==space),max(v['coverage'] for v in intervals if v['space']==space)] for space in ['response','treated']},intervalMeanWidthRange={space:[min(v['meanWidth'] for v in intervals if v['space']==space),max(v['meanWidth'] for v in intervals if v['space']==space)] for space in ['response','treated']},maximumPointError=max(c['maximumPointError'] for c in comparisons),maximumIntervalBoundError=max(v['maximumBoundError'] for v in intervals),modelChecks=models,inputFreezeSHA256=sha(src/'input-freeze.json'),predictionFreezeSHA256=sha(native/'prediction-freeze.json'),scorerSHA256=sha(Path(__file__)),packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','anndata']},scope='All held-out population RNA donor-times in an already inspected single study; no independent external validation, general uncertainty calibration, unseen perturbation, phenotype or clinical qualification.')
 write(a.out/'summary.json',summary);write(a.out/'scores.json',scores);write(a.out/'intervals.json',intervals);write(a.out/'comparisons.json',comparisons);np.savez_compressed(a.out/'all-vectors.npz',featureIDs=np.array(panel),**arrays);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
