#!/usr/bin/env python3
"""Independent complete-panel reconstruction and scoring of frozen native outputs."""
import argparse,hashlib,json,importlib.metadata
from pathlib import Path
import numpy as np
from sklearn.kernel_ridge import KernelRidge

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--native',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 input_freeze=json.loads((a.inputs/'input-freeze.json').read_text());prediction_freeze=json.loads((a.native/'prediction-freeze.json').read_text())
 for name,digest in input_freeze['files'].items():assert sha(a.inputs/name)==digest,name
 # Full native trees may be retained remotely. Every scientific JSON used here
 # must be present and match its pre-scoring freeze; absent H5AD copies are not reverified here.
 for name,digest in prediction_freeze['files'].items():
  if (a.native/name).is_file():assert sha(a.native/name)==digest,name
 assert prediction_freeze['completedFolds']==26 and prediction_freeze['scoringStarted'] is False
 execution=json.loads((a.native/'execution.json').read_text())
 assert execution['inputFreezeSHA256']==sha(a.inputs/'input-freeze.json')
 panel=json.loads((a.inputs/'panel.json').read_text());cohorts={s:json.loads((a.inputs/(s+'-cohort.json')).read_text()) for s in ['Kang','HIRISA']};counts={};logs={};ids={}
 for s in cohorts:
  with np.load(a.inputs/(s+'-source-counts.npz'),allow_pickle=False) as f:counts[s]=f['counts'].copy();ids[s]=f['featureIDs'].tolist()
  logs[s]=np.log1p(counts[s].astype(np.float64)/counts[s].sum(axis=1,keepdims=True,dtype=np.uint64)*1e6)
 order={s:np.array([ids[s].index(x) for x in panel]) for s in cohorts};checks=[];scores=[]
 for fold in json.loads((a.inputs/'folds.json').read_text()):
  for suffix in ['model/model.json','model/training/report.json','prediction/report.json']:
   name=fold['id']+'/'+suffix
   assert name in prediction_freeze['files'] and sha(a.native/name)==prediction_freeze['files'][name],name
  root=a.native/fold['id'];model=json.loads((root/'model/model.json').read_text());report=json.loads((root/'prediction/report.json').read_text());prediction=report['predictions'][0]
  assert model['featureIDs']==report['featureIDs']==panel and len(report['predictions'])==1
  train=fold['trainingStudy'];query=fold['queryStudy'];samples=cohorts[train]['samples'];train_indices=fold['trainingIndices'];donors=sorted({samples[i]['donorID'] for i in train_indices});assert donors==model['trainingDonors'] and fold['heldOutDonor'] not in donors
  native_bulk=json.loads((root/'model/training/report.json').read_text())['pseudobulk']
  assert native_bulk['featureIDs']==ids[train] and len(native_bulk['groups'])==len(train_indices)
  for row,group in enumerate(native_bulk['groups']):
   matches=[i for i in train_indices if group['sampleIDs']==[samples[i]['id']]];assert len(matches)==1
   m=native_bulk['matrix'];begin,end=m['rowOffsets'][row:row+2];values=np.zeros(len(ids[train]),dtype=np.uint64)
   values[m['featureIndices'][begin:end]]=m['counts'][begin:end]
   np.testing.assert_array_equal(values,counts[train][matches[0]])
  pairs=[]
  for donor in donors:
   c=[i for i in train_indices if samples[i]['donorID']==donor and samples[i]['condition']=='control'];t=[i for i in train_indices if samples[i]['donorID']==donor and samples[i]['condition']=='IFNB'];assert len(c)==len(t)==1;pairs.append((c[0],t[0]))
  c,t=np.array(pairs).T;controls=logs[train][c][:,order[train]];response=logs[train][t][:,order[train]]-controls
  pc=counts[train][c][:,order[train]]+counts[train][t][:,order[train]]
  selected=np.flatnonzero((pc.sum(axis=0,dtype=np.uint64)>=10)&((pc>0).sum(axis=0)>=2));assert selected.tolist()==model['selectedFeatureIndices']
  center=controls[:,selected].mean(axis=0);scale=controls[:,selected].std(axis=0);constant=np.all(controls[:,selected]==controls[0,selected],axis=0);center[constant]=controls[0,selected][constant];scale[constant|(scale==0)]=1
  context=(controls[:,selected]-center)/scale/np.sqrt(len(selected));mean=response.mean(axis=0);median=np.median(response,axis=0)
  kernel=context@context.T;dual=np.linalg.solve(kernel+np.eye(len(donors)),response-mean)
  estimator=KernelRidge(alpha=1,kernel='precomputed').fit(kernel,response-mean)
  np.testing.assert_allclose(estimator.dual_coef_,dual,rtol=1e-8,atol=1e-9)
  for key,expected in [('contextCenters',center),('contextScales',scale),('contexts',context),('meanResponse',mean),('medianResponse',median),('dualCoefficients',dual)]:np.testing.assert_allclose(model[key],expected,rtol=1e-8,atol=1e-9,err_msg=fold['id']+'/'+key)
  control=logs[query][fold['queryIndex'],order[query]];np.testing.assert_allclose(prediction['control'],control,rtol=1e-10,atol=1e-10)
  assert prediction['libraryCounts']==int(counts[query][fold['queryIndex']].sum(dtype=np.uint64))
  qc=(control[selected]-center)/scale/np.sqrt(len(selected));ridge=estimator.predict((qc@context.T)[None,:])[0]+mean
  refs=dict(noChange=np.zeros(len(panel)),meanResponse=mean,medianResponse=median,contextRidge=ridge)
  assert [x['baseline'] for x in prediction['estimates']]==list(refs)
  truth=logs[query][fold['scoringIndex'],order[query]]
  # Use the same native control for observed and predicted changes so no-change
  # stays exactly zero rather than reflecting cross-library roundoff.
  native_control=np.array(prediction['control']);observed=truth-native_control
  for e in prediction['estimates']:
   baseline=e['baseline'];expected=np.maximum(0,control+refs[baseline]);actual=np.array(e['predictedTreated']);change=actual-native_control
   for key,value in [('unclippedResponse',refs[baseline]),('predictedTreated',expected),('predictedResponse',expected-control)]:np.testing.assert_allclose(e[key],value,rtol=1e-8,atol=1e-8,err_msg=fold['id']+'/'+baseline+'/'+key)
   np.testing.assert_allclose(e['impliedCPMSum'],np.expm1(expected).sum(),rtol=1e-9,atol=1e-5)
   error=actual-truth;correlation=None if np.std(change)==0 or np.std(observed)==0 else float(np.corrcoef(change,observed)[0,1])
   scores.append(dict(fold=fold['id'],mode=fold['mode'],trainingStudy=train,queryStudy=query,donor=fold['heldOutDonor'],baseline=baseline,features=len(panel),responseRMSE=float(np.sqrt(np.mean(error**2))),responseMAE=float(np.mean(np.abs(error))),responsePearson=correlation,panelImpliedCPM=e['impliedCPMSum']))
   checks.append(dict(fold=fold['id'],baseline=baseline,maximumPredictionDifference=float(np.max(np.abs(actual-expected))),maximumDualDifference=float(np.max(np.abs(np.array(model['dualCoefficients'])-dual))),contextFeatures=len(selected)))
 summaries=[]
 for target in ['Kang','HIRISA']:
  subset=[s for s in scores if s['queryStudy']==target]
  means={mode:{b:float(np.mean([s['responseRMSE'] for s in subset if s['mode']==mode and s['baseline']==b])) for b in refs} for mode in ['cross','within']}
  donor_scores={d:{s['baseline']:s['responseRMSE'] for s in subset if s['mode']=='cross' and s['donor']==d} for d in sorted({s['donor'] for s in subset})}
  cross=means['cross'];summaries.append(dict(queryStudy=target,donors=len(donor_scores),means=means,primaryPass=cross['contextRidge']<min(cross['noChange'],cross['meanResponse']),ridgeGainOverMeanPercent=100*(1-cross['contextRidge']/cross['meanResponse']),ridgeWorseNoChange=sum(x['contextRidge']>x['noChange'] for x in donor_scores.values()),ridgeWorseMean=sum(x['contextRidge']>x['meanResponse'] for x in donor_scores.values())))
 write(a.out/'scores.json',scores);write(a.out/'comparisons.json',checks);write(a.out/'summary.json',summaries)
 write(a.out/'checks.json',dict(status='passed',folds=26,vectors=len(checks),features=len(panel),maximumPredictionDifference=max(x['maximumPredictionDifference'] for x in checks),maximumDualDifference=max(x['maximumDualDifference'] for x in checks),inputFreezeSHA256=sha(a.inputs/'input-freeze.json'),predictionFreezeSHA256=sha(a.native/'prediction-freeze.json'),scorerSHA256=sha(__file__),packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','scikit-learn']}))
 print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
