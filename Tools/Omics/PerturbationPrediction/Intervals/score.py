#!/usr/bin/env python3
"""Independently check future-donor intervals, point preservation and empirical coverage."""
import argparse,hashlib,json,tarfile,importlib.metadata
from pathlib import Path
import numpy as np
from scipy.stats import t

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,sort_keys=True,indent=2,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--native',type=Path,required=True);p.add_argument('--previous-native',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 prior=json.loads((a.inputs/'input-freeze.json').read_text());input_freeze=json.loads((a.native/'input-freeze.json').read_text());freeze=json.loads((a.native/'prediction-freeze.json').read_text());previous_freeze=json.loads((a.previous_native/'prediction-freeze.json').read_text())
 assert input_freeze['priorInputFreezeSHA256']==sha(a.inputs/'input-freeze.json') and freeze['inputFreezeSHA256']==sha(a.native/'input-freeze.json') and not freeze['scoringStarted']
 for name,digest in prior['files'].items():assert sha(a.inputs/name)==digest,name
 for name,digest in input_freeze['files'].items():assert sha(a.native/name)==digest,name
 folds=json.loads((a.inputs/'folds.json').read_text());by_id={f['id']:f for f in folds};panel=json.loads((a.inputs/'panel.json').read_text());counts={};logs={};order={};cohorts={}
 for s in ['Kang','HIRISA']:
  with np.load(a.inputs/(s+'-source-counts.npz'),allow_pickle=False) as f:counts[s]=f['counts'].copy();ids=f['featureIDs'].tolist()
  logs[s]=np.log1p(counts[s].astype(np.float64)/counts[s].sum(axis=1,keepdims=True,dtype=np.uint64)*1e6);order[s]=[ids.index(x) for x in panel];cohorts[s]=json.loads((a.inputs/(s+'-cohort.json')).read_text())
 records=[];checks=[]
 for archive in freeze['folds']:
  fold=by_id[archive['fold']];path=a.native/archive['path'];assert path.stat().st_size==archive['bytes'] and sha(path)==archive['SHA256']
  with tarfile.open(path,'r:gz') as bundle:
   members=bundle.getmembers();assert len(members)==len(archive['members']) and {x.name for x in members}==set(archive['members'])
   for member in members:
    assert member.isfile();b=bundle.extractfile(member).read();expected=archive['members'][member.name];assert len(b)==expected['bytes'] and hashlib.sha256(b).hexdigest()==expected['SHA256']
   model=json.load(bundle.extractfile('model/model.json'));prediction=json.load(bundle.extractfile('prediction/report.json'))['predictions'][0]
  old_root=a.previous_native/fold['id']
  for suffix in ['model/model.json','prediction/report.json']:
   key=fold['id']+'/'+suffix;assert sha(old_root/suffix)==previous_freeze['files'][key]
  old_model=json.loads((old_root/'model/model.json').read_text());old_prediction=json.loads((old_root/'prediction/report.json').read_text())['predictions'][0]
  assert prediction['control']==old_prediction['control'] and prediction['estimates']==old_prediction['estimates']
  for key in ['featureIDs','meanResponse','medianResponse','selectedFeatureIndices','contextCenters','contextScales','contexts','dualCoefficients','maximumSolveResidual']:assert model[key]==old_model[key],key
  source=fold['trainingStudy'];target=fold['queryStudy'];samples=cohorts[source]['samples'];donors=model['trainingDonors'];rows=fold['trainingIndices'];pairs=[]
  for donor in donors:
   c=[i for i in rows if samples[i]['donorID']==donor and samples[i]['condition']=='control'];r=[i for i in rows if samples[i]['donorID']==donor and samples[i]['condition']=='IFNB'];assert len(c)==len(r)==1;pairs.append((c[0],r[0]))
  c,r=np.array(pairs).T;responses=logs[source][r][:,order[source]]-logs[source][c][:,order[source]];mean=responses.mean(axis=0);variance=responses.var(axis=0,ddof=1);available=(~np.all(responses==responses[0],axis=0))&(variance>0)
  native_variance=np.array([np.nan if v is None else v for v in model['donorResponseVariances']]);assert np.array_equal(np.isfinite(native_variance),available)
  np.testing.assert_allclose(native_variance[available],variance[available],atol=1e-11,rtol=1e-9)
  interval=prediction['meanResponsePredictiveInterval'];n=len(donors);critical=float(t.ppf(.975,n-1));assert interval['nominalCoverage']==.95 and interval['degreesOfFreedom']==n-1 and interval['trainingDonors']==n
  np.testing.assert_allclose(interval['studentCriticalValue'],critical,atol=1e-11,rtol=1e-11)
  assert interval['unavailableFeatureIndices']==np.flatnonzero(~available).tolist()
  control=np.array(prediction['control']);half=critical*np.sqrt(variance*(1+1/n));lower=mean-half;upper=mean+half
  refs={'unclippedResponseLower':lower,'unclippedResponseUpper':upper,'predictedTreatedLower':np.maximum(0,control+lower),'predictedTreatedUpper':np.maximum(0,control+upper)};bounds={};maximum=0.
  for key,ref in refs.items():
   values=np.array([np.nan if x is None else x for x in interval[key]]);assert np.array_equal(np.isfinite(values),available);np.testing.assert_allclose(values[available],ref[available],atol=1e-9,rtol=1e-9);bounds[key]=values;maximum=max(maximum,float(np.max(np.abs(values[available]-ref[available]))) if available.any() else 0.)
  truth=logs[target][fold['scoringIndex'],order[target]]
  for space,observation,lo,hi in [('unclippedResponse',truth-control,bounds['unclippedResponseLower'],bounds['unclippedResponseUpper']),('predictedTreated',truth,bounds['predictedTreatedLower'],bounds['predictedTreatedUpper'])]:
   y=observation[available];l=lo[available];h=hi[available];covered=(y>=l)&(y<=h)
   records.append(dict(fold=fold['id'],mode=fold['mode'],queryStudy=target,donor=fold['heldOutDonor'],space=space,nominalCoverage=.95,trainingDonors=n,totalFeatures=len(panel),availableFeatures=int(available.sum()),unavailableFeatures=int((~available).sum()),coveredFeatures=int(covered.sum()),belowFeatures=int((y<l).sum()),aboveFeatures=int((y>h).sum()),empiricalCoverage=float(covered.mean()) if len(y) else None,meanWidth=float(np.mean(h-l)) if len(y) else None))
  checks.append(dict(fold=fold['id'],availableFeatures=int(available.sum()),maximumBoundDifference=maximum,maximumVarianceDifference=float(np.max(np.abs(native_variance[available]-variance[available]))) if available.any() else 0.,studentCriticalDifference=abs(interval['studentCriticalValue']-critical),fourPointPredictionsExact=True))
 summaries=[]
 for target in ['Kang','HIRISA']:
  for mode in ['within','cross']:
   for space in ['unclippedResponse','predictedTreated']:
    selected=[r for r in records if r['queryStudy']==target and r['mode']==mode and r['space']==space];assert len(selected)==(8 if target=='Kang' else 5)
    summaries.append(dict(queryStudy=target,mode=mode,space=space,donors=len(selected),nominalCoverage=.95,meanDonorCoverage=float(np.mean([r['empiricalCoverage'] for r in selected])),minimumDonorCoverage=min(r['empiricalCoverage'] for r in selected),maximumDonorCoverage=max(r['empiricalCoverage'] for r in selected),meanDonorWidth=float(np.mean([r['meanWidth'] for r in selected])),meanAvailableFeatureFraction=float(np.mean([r['availableFeatures']/r['totalFeatures'] for r in selected]))))
 write(a.out/'scores.json',records);write(a.out/'comparisons.json',checks);write(a.out/'summary.json',summaries)
 write(a.out/'checks.json',dict(status='passed',folds=len(checks),intervalBoundArrays=len(checks)*4,unchangedPointVectors=len(checks)*4,maximumBoundDifference=max(x['maximumBoundDifference'] for x in checks),maximumVarianceDifference=max(x['maximumVarianceDifference'] for x in checks),maximumStudentCriticalDifference=max(x['studentCriticalDifference'] for x in checks),scorerSHA256=sha(__file__),inputFreezeSHA256=sha(a.native/'input-freeze.json'),predictionFreezeSHA256=sha(a.native/'prediction-freeze.json'),packages={n:importlib.metadata.version(n) for n in ['numpy','scipy']}))
 print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
