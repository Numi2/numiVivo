"""Independently reconstruct every full and donor-omission paired result.

Uses the original NumPy cell moments, not native fitted calibration parameters.
LAPACK symmetric eigenvalues provide a separate covariance calculation.
"""
from pathlib import Path
import collections,gzip,hashlib,json,sys
import numpy as np

root=Path(sys.argv[1]);origin=sys.argv[2]
parent=Path(sys.argv[3]) if len(sys.argv)>3 else Path('/Users/home/numivivo-count-calibration-20260911')
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save_gzip(path,obj):
 with path.open('wb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode())
freeze=read(root/'input-freeze.json')['origins'][origin]
input_bytes=gzip.decompress((root/(origin+'-input.json.gz')).read_bytes());inp=json.loads(input_bytes)
assert hashlib.sha256(input_bytes).hexdigest()==freeze['inputSHA256']
ref_path=parent/(origin+'-reference-moments.json.gz');reference=read(ref_path)
meta=read(parent/'inputs'/(origin+'.json'));donors=sorted({g['donorID'] for g in meta['groups']});features=inp['featureIDs'];m=len(features)
assert sha(parent/'inputs'/(origin+'.json'))==freeze['parentMetadataSHA256']
assert len(reference)==len(meta['groups'])==len(inp['groups'])
for original,g,r in zip(meta['groups'],inp['groups'],reference):
 assert original['donorID']==g['donorID'] and original['conditionID']==g['conditionID']
 for k in ['cells','libraryCounts','counts','positiveCells']:assert g['moments'][k]==r[k]
errors=[];maximum_error=0.;values_checked=0;expected_reports=[];all_native=[];fold_summaries=[]
def marginal(groups):
 n=len(groups);N=np.array([g['cells'] for g in groups],float)[:,None]
 means=np.array([g['meanCPM'] for g in groups]);variances=np.array([g['sampleVarianceCPM'] for g in groups]);shot=np.array([g['meanPoissonVarianceCPM'] for g in groups]);pairs=np.array([g['distinctCellRateProductCPM2'] for g in groups])
 denominator=((N-1)*pairs).sum(0);numerator=((N-1)*(variances-shot)).sum(0)
 phi=np.maximum(0,np.divide(numerator,denominator,out=np.zeros(m),where=denominator>0))
 admitted=(denominator>0)&((phi==0)|((phi>=1e-8)&(phi<=100)))
 mean=means.mean(0);observed=means.var(0,ddof=1);measurement=((shot+phi*pairs)/N).mean(0)
 future=np.maximum(0,observed-measurement)*(1+1/n)+measurement/n
 shape=np.divide(mean*mean,future,out=np.zeros(m),where=future>0);rate=np.divide(mean,future,out=np.zeros(m),where=future>0)
 statuses=np.full(m,'insufficientWithinDonorCountPairs',dtype=object)
 statuses[(denominator>0)&~admitted]='cellDispersionOutsideObservationKernelDomain'
 statuses[admitted]='unidentifiedPositiveRatePrior'
 prior=admitted&(mean>0)&(future>0)
 statuses[prior]='gammaPriorOutsideObservationKernelDomain'
 statuses[prior&(shape>=.1)&(shape<=1e6)&(rate>=1e-12)&(rate<=1e12)]='availableConditionalMomentCalibration'
 return means,mean,observed,measurement,admitted,statuses
with gzip.open(root/(origin+'-native.jsonl.gz'),'rt') as f:
 for row_number,line in enumerate(f):
  obj=json.loads(line);report=obj['report'];excluded=obj.get('excludedDonorID');expected_excluded=([None]+donors)[row_number]
  assert obj['origin']==origin and excluded==expected_excluded
  assert bytes(report['trainingSource']['bytes']).hex()==freeze['inputSHA256']
  ids=[d for d in donors if d!=excluded];assert report['donorIDs']==ids
  assert report['controlConditionID']=='control' and report['treatedConditionID']=='IFNB'
  assert report['method']=='paired-donor-count-moment-compatibility-v1' and report['covariancePSDRelativeTolerance']==1e-10
  groups=[[reference[next(i for i,g in enumerate(meta['groups']) if g['donorID']==d and g['conditionID']==condition)] for d in ids] for condition in ['control','IFNB']]
  for prefix,gs in zip(['control','treated'],groups):
   assert report[prefix+'Cells']==[g['cells'] for g in gs]
   assert report[prefix+'LibraryCounts']==[g['libraryCounts'] for g in gs]
  c,mc,vc,nc,ac,sc=marginal(groups[0]);t,mt,vt,nt,at,st=marginal(groups[1]);n=len(ids)
  covariance=((c-mc)*(t-mt)).sum(0)/(n-1);response=(t-c).var(0,ddof=1)
  bulk=[np.array([g['counts'] for g in gs])/np.array([g['libraryCounts'] for g in gs])[:,None]*1e6 for gs in groups]
  bc,bt=bulk;valid=ac&at;rc=vc-nc;rt=vt-nt
  # Batched independent symmetric eigensolver, retaining signed eigenvalues.
  raw=np.zeros((m,2,2));raw[:,0,0]=rc;raw[:,1,1]=rt;raw[:,0,1]=raw[:,1,0]=covariance
  clipped=raw.copy();clipped[:,0,0]=np.maximum(0,rc);clipped[:,1,1]=np.maximum(0,rt)
  re=np.linalg.eigvalsh(raw)[:,0];ce=np.linalg.eigvalsh(clipped)[:,0]
  rs=np.where(re < -1e-10*np.abs(raw).max((1,2)),'indefinite','positiveSemidefinite').astype(object)
  cs=np.where(ce < -1e-10*np.abs(clipped).max((1,2)),'indefinite','positiveSemidefinite').astype(object)
  rs[~valid]=cs[~valid]='unavailableCellDispersion'
  expected={'controlMeanCPM':mc,'treatedMeanCPM':mt,'observedControlVarianceCPM2':vc,'observedTreatedVarianceCPM2':vt,'observedCovarianceCPM2':covariance,'observedResponseVarianceCPM2':response,
   'controlMeasurementVarianceCPM2':np.where(ac,nc,np.nan),'treatedMeasurementVarianceCPM2':np.where(at,nt,np.nan),
   'rawLatentControlVarianceCPM2':np.where(valid,rc,np.nan),'rawLatentTreatedVarianceCPM2':np.where(valid,rt,np.nan),'rawLatentResponseVarianceCPM2':np.where(valid,response-nc-nt,np.nan),
   'rawMinimumEigenvalueCPM2':np.where(valid,re,np.nan),'marginalClippedMinimumEigenvalueCPM2':np.where(valid,ce,np.nan),
   'meanCellRateResponseCPM':mt-mc,'meanPseudobulkResponseCPM':(bt-bc).mean(0),'meanLog1pCellRateResponse':(np.log1p(t)-np.log1p(c)).mean(0),'meanLog1pPseudobulkResponse':(np.log1p(bt)-np.log1p(bc)).mean(0),
   'maximumAbsoluteEndpointDifferenceCPM':np.maximum(np.abs(c-bc).max(0),np.abs(t-bt).max(0)),
   'maximumAbsoluteLog1pEndpointDifference':np.maximum(np.abs(np.log1p(c)-np.log1p(bc)).max(0),np.abs(np.log1p(t)-np.log1p(bt)).max(0))}
  fs=report['features'];assert [f['featureID'] for f in fs]==features
  for name,array in expected.items():
   actual=np.array([f.get(name,np.nan) for f in fs]);mask=np.isfinite(array)
   assert np.array_equal(np.isfinite(actual),mask),(excluded,name,'missingness')
   err=float(np.max(np.abs(actual[mask]-array[mask])/np.maximum(1,np.abs(array[mask])))) if mask.any() else 0.
   maximum_error=max(maximum_error,err);values_checked+=int(mask.sum())
   if err>2e-7 or not np.isfinite(err):errors.append(dict(excluded=excluded,field=name,maximumScaledError=err))
  for name,array in [('controlCalibrationStatus',sc),('treatedCalibrationStatus',st),('rawCovarianceStatus',rs),('marginalClippedCovarianceStatus',cs)]:
   assert [f[name] for f in fs]==array.tolist(),(excluded,name)
  endpoint=expected['meanLog1pCellRateResponse']-expected['meanLog1pPseudobulkResponse']
  summary={'excludedDonorID':excluded,'donors':n,'features':m,'rawCovarianceStates':dict(collections.Counter(rs)),'marginalClippedCovarianceStates':dict(collections.Counter(cs)),
   'bothGammaCalibrationsAvailable':int(((sc=='availableConditionalMomentCalibration')&(st=='availableConditionalMomentCalibration')).sum()),
   'negativeObservedCovariance':int((covariance<0).sum()),'negativeCorrectedResponseVariance':int((valid&(response-nc-nt<0)).sum()),
   'logEndpointResponseSignReversals':int((expected['meanLog1pCellRateResponse']*expected['meanLog1pPseudobulkResponse']<0).sum()),
   'medianAbsoluteLogResponseEndpointDifference':float(np.median(np.abs(endpoint))),'maximumAbsoluteLogResponseEndpointDifference':float(np.abs(endpoint).max()),
   'maximumDonorConditionLogEndpointDifference':float(expected['maximumAbsoluteLog1pEndpointDifference'].max())}
  fold_summaries.append(summary)
  expected_reports.append({'excludedDonorID':excluded,'values':{k:[float(x) if np.isfinite(x) else None for x in a] for k,a in expected.items()},'rawCovarianceStatus':rs.tolist(),'marginalClippedCovarianceStatus':cs.tolist()})
  all_native.append(fs)
assert len(all_native)==len(donors)+1
sens=[]
for j,feature in enumerate(features):
 fs=[fold[j] for fold in all_native];cov=[f['observedCovarianceCPM2'] for f in fs];raw=[f['rawCovarianceStatus'] for f in fs];clipped=[f['marginalClippedCovarianceStatus'] for f in fs]
 sens.append({'featureID':feature,'covarianceChangesSignAcrossDonorOmissions':min(cov)<0<max(cov),'rawStates':dict(collections.Counter(raw)),'marginalClippedStates':dict(collections.Counter(clipped)),
  'fullRawState':raw[0],'fullClippedState':clipped[0],'maximumAbsoluteOmissionLogResponseChange':max(abs(f['meanLog1pCellRateResponse']-fs[0]['meanLog1pCellRateResponse']) for f in fs[1:])})
save_gzip(root/(origin+'-reference-results.json.gz'),expected_reports);save_gzip(root/(origin+'-sensitivity.json.gz'),sens)
result={'status':'passed' if not errors else 'failed','origin':origin,'reports':len(all_native),'featureRecords':m*len(all_native),'numericalValuesCompared':values_checked,'maximumScaledError':maximum_error,'errors':errors,
 'independentReferenceMomentsSHA256':sha(ref_path),'inputSHA256':freeze['inputSHA256'],'folds':fold_summaries,
 'sensitivity':{'covarianceSignChanges':sum(s['covarianceChangesSignAcrossDonorOmissions'] for s in sens),'rawStateChanges':sum(len(s['rawStates'])>1 for s in sens),'clippedStateChanges':sum(len(s['marginalClippedStates'])>1 for s in sens)},
 'interpretation':'Training-donor sensitivity, not leave-one-donor-out prediction or independent biological validation.'}
(root/(origin+'-verification.json')).write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
print(json.dumps(result));sys.exit(bool(errors))
