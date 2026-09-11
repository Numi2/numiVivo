"""Independent NumPy calibration from original per-cell moments; no native parameters reused."""
from pathlib import Path
import collections,gzip,hashlib,json,sys
import numpy as np
root=Path(sys.argv[1]);origin=sys.argv[2]
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
meta=read(root/'inputs'/(origin+'.json'));native=read(root/(origin+'-native.json.gz'));reference=read(root/(origin+'-reference-moments.json.gz'))
errors=[];maxmoment=0.;maxcal=0.;status=[]
def close(a,b,tol,where):
 global maxmoment,maxcal
 a=np.asarray(a,float);b=np.asarray(b,float);assert a.shape==b.shape,where
 e=float(np.max(np.abs(a-b)/np.maximum(1,np.abs(b)))) if a.size else 0
 if where.startswith('moment'):maxmoment=max(maxmoment,e)
 else:maxcal=max(maxcal,e)
 if not np.isfinite(e) or e>tol:errors.append({'where':where,'maximumScaledError':e})
assert len(native['moments'])==len(reference)==len(meta['groups'])
for i,(a,b) in enumerate(zip(native['moments'],reference)):
 for key in ['cells','libraryCounts','counts','positiveCells']:assert a[key]==b[key],(i,key)
 for key in ['meanCPM','sampleVarianceCPM','meanPoissonVarianceCPM','distinctCellRateProductCPM2']:close(a[key],b[key],2e-8,f'moment/{i}/{key}')
for model in native['models']:
 condition=model['conditionID'];indices=[i for i,g in enumerate(meta['groups']) if g['conditionID']==condition]
 assert model['trainingDonorIDs']==[meta['groups'][i]['donorID'] for i in indices]
 assert bytes(model['trainingSource']['bytes']).hex()==native['bindingSHA256']
 assert [f['featureID'] for f in model['features']]==meta['featureIDs']
 groups=[reference[i] for i in indices];N=np.array([g['cells'] for g in groups],float)[:,None];d=len(groups)
 means=np.array([g['meanCPM'] for g in groups]);variances=np.array([g['sampleVarianceCPM'] for g in groups]);shot=np.array([g['meanPoissonVarianceCPM'] for g in groups]);pairs=np.array([g['distinctCellRateProductCPM2'] for g in groups])
 numer=((N-1)*(variances-shot)).sum(axis=0);denom=((N-1)*pairs).sum(axis=0);mean=means.mean(axis=0);observed=means.var(axis=0,ddof=1)
 for j,f in enumerate(model['features']):
  expected=dict(cellDispersionNumerator=numer[j],cellDispersionDenominator=denom[j],meanDonorRateCPM=mean[j],observedDonorRateVarianceCPM=observed[j])
  state='insufficientWithinDonorCountPairs';poisson=False;boundary=False
  if denom[j]>0:
   raw=numer[j]/denom[j];phi=max(0,raw);expected['rawCellDispersion']=raw
   if np.isfinite(phi) and (phi==0 or 1e-8<=phi<=100):
    mv=float(((shot[:,j]+phi*pairs[:,j])/N[:,0]).mean());latent=max(0,observed[j]-mv);future=latent*(1+1/d)+mv/d
    expected.update(cellDispersion=phi,meanDonorMeasurementVarianceCPM=mv,latentDonorRateVarianceCPM=latent,newDonorRateVarianceCPM=future)
    poisson=raw<=0;boundary=observed[j]-mv<=0
    if mean[j]>0 and future>0:
     a=mean[j]**2/future;b=mean[j]/future
     if np.isfinite(a) and np.isfinite(b) and .1<=a<=1e6 and 1e-12<=b<=1e12:
      expected.update(gammaPriorShape=a,gammaPriorRatePerCPM=b);state='availableConditionalMomentCalibration'
     else:state='gammaPriorOutsideObservationKernelDomain'
    else:state='unidentifiedPositiveRatePrior'
   else:state='cellDispersionOutsideObservationKernelDomain'
  assert f['status']==state,(j,f['status'],state)
  assert f['poissonBoundary']==bool(poisson) and f['latentVarianceBoundary']==bool(boundary)
  for key,value in expected.items():close(f[key],value,2e-7,f'{condition}/{j}/{key}')
  optional=['rawCellDispersion','cellDispersion','meanDonorMeasurementVarianceCPM','latentDonorRateVarianceCPM','newDonorRateVarianceCPM','gammaPriorShape','gammaPriorRatePerCPM']
  assert all(k in f if k in expected else k not in f for k in optional)
  assert f['trainingCounts']==sum(g['counts'][j] for g in groups)
  assert f['positiveDonors']==sum(g['counts'][j]>0 for g in groups)
 rates=np.array([g['depthRateCorrelation'] for g in groups]);support=np.array([g['positiveCells'] for g in groups])>=10
 status.append({'condition':condition,'features':len(model['features']),'states':dict(collections.Counter(f['status'] for f in model['features'])),'poissonBoundary':sum(f['poissonBoundary'] for f in model['features']),'latentVarianceBoundary':sum(f['latentVarianceBoundary'] for f in model['features']),'depthRateAssociation':{'supportedDonorGenes':int(support.sum()),'absolutePearsonAbovePoint2':int(((np.abs(rates)>.2)&support).sum()),'medianAbsolutePearson':float(np.median(np.abs(rates[support]))),'interpretation':'Descriptive RNA-depth association; mathematical coupling and biological mixtures are not separated. No association-based feature exclusion.'}})
result={'status':'passed' if not errors else 'failed','origin':origin,'sourceCells':native['cells'],'sourceNonzeros':native['nonzeros'],'groups':len(reference),'fullRNAFeatures':len(meta['featureIDs']),'momentValuesCompared':len(reference)*len(meta['featureIDs'])*4,'maximumScaledMomentError':maxmoment,'maximumScaledCalibrationError':maxcal,'comparisonScale':'abs(native-reference)/max(1,abs(reference))','sourceBinding':native['bindingSHA256'],'models':status,'errors':errors}
(root/(origin+'-verification.json')).write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
print(json.dumps(result));sys.exit(bool(errors))
