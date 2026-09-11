"""Independent NumPy reconstruction of every training-only calibration field."""
from pathlib import Path
import collections,gzip,hashlib,json,sys,time
import numpy as np
root=Path(sys.argv[1]);full=Path(sys.argv[2]);protocol=json.loads((root/'protocol.json').read_text());assert json.loads((root/'state.json').read_text())['status']=='completed-all-training-only-folds'
values=0;maximum=0.;folds=[];started=time.time()
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def close(a,b,label):
 global values,maximum
 e=abs(float(a)-float(b))/max(1,abs(float(b)));assert np.isfinite(e) and e<=2e-7,(label,a,b,e);values+=1;maximum=max(maximum,e)
for fold in protocol['folds']:
 tag=fold['tag'];compressed=(root/'inputs'/(tag+'.json.gz')).read_bytes();raw=gzip.decompress(compressed);inp=json.loads(raw);parent=read(root/(fold['origin']+'-parent-input.json.gz'));parent_raw=gzip.decompress((root/(fold['origin']+'-parent-input.json.gz')).read_bytes())
 assert hashlib.sha256(raw).hexdigest()==fold['inputSHA256'] and hashlib.sha256(compressed).hexdigest()==fold['compressedInputSHA256'];assert hashlib.sha256(parent_raw).hexdigest()==inp['parentInputSHA256']==fold['parentInputSHA256']
 selected=[g for g in parent['groups'] if g['donorID']!=fold['excludedDonorID']];assert inp['groups']==selected and inp['featureIDs']==parent['featureIDs'];assert not any(g['donorID']==fold['excludedDonorID'] for g in inp['groups'])
 receipt=read(root/'native'/(tag+'-receipt.json'));compressed=(root/'native'/(tag+'.json.gz')).read_bytes();raw=gzip.decompress(compressed);native=json.loads(raw);assert hashlib.sha256(compressed).hexdigest()==receipt['compressedOutputSHA256'] and hashlib.sha256(raw).hexdigest()==receipt['outputSHA256'];assert native['inputSHA256']==fold['inputSHA256']
 models={};states={}
 for side,condition in [('control',inp['controlConditionID']),('treated',inp['treatedConditionID'])]:
  model=native[side];ids=sorted({g['donorID'] for g in selected});assert ids==model['trainingDonorIDs']==fold['trainingDonorIDs'];assert bytes(model['trainingSource']['bytes']).hex()==fold['inputSHA256'];assert [f['featureID'] for f in model['features']]==inp['featureIDs'] and model['conditionID']==condition
  groups=[next(g['moments'] for g in selected if g['donorID']==d and g['conditionID']==condition) for d in ids];assert model['trainingCellCounts']==[g['cells'] for g in groups]
  N=np.array([g['cells'] for g in groups],float)[:,None];d=len(groups);means=np.array([g['meanCPM'] for g in groups]);variances=np.array([g['sampleVarianceCPM'] for g in groups]);shot=np.array([g['meanPoissonVarianceCPM'] for g in groups]);pairs=np.array([g['distinctCellRateProductCPM2'] for g in groups]);numer=((N-1)*(variances-shot)).sum(axis=0);denom=((N-1)*pairs).sum(axis=0);mean=means.mean(axis=0);observed=means.var(axis=0,ddof=1)
  for j,f in enumerate(model['features']):
   expected=dict(cellDispersionNumerator=numer[j],cellDispersionDenominator=denom[j],meanDonorRateCPM=mean[j],observedDonorRateVarianceCPM=observed[j]);state='insufficientWithinDonorCountPairs';poisson=False;boundary=False
   if denom[j]>0:
    rawphi=numer[j]/denom[j];phi=max(0,rawphi);expected['rawCellDispersion']=rawphi
    if np.isfinite(phi) and (phi==0 or 1e-8<=phi<=100):
     mv=float(((shot[:,j]+phi*pairs[:,j])/N[:,0]).mean());latent=max(0,observed[j]-mv);future=latent*(1+1/d)+mv/d;expected.update(cellDispersion=phi,meanDonorMeasurementVarianceCPM=mv,latentDonorRateVarianceCPM=latent,newDonorRateVarianceCPM=future);poisson=rawphi<=0;boundary=observed[j]-mv<=0
     if mean[j]>0 and future>0:
      a=mean[j]**2/future;b=mean[j]/future
      if np.isfinite(a) and np.isfinite(b) and .1<=a<=1e6 and 1e-12<=b<=1e12:expected.update(gammaPriorShape=a,gammaPriorRatePerCPM=b);state='availableConditionalMomentCalibration'
      else:state='gammaPriorOutsideObservationKernelDomain'
     else:state='unidentifiedPositiveRatePrior'
    else:state='cellDispersionOutsideObservationKernelDomain'
   assert f['status']==state and f['poissonBoundary']==bool(poisson) and f['latentVarianceBoundary']==bool(boundary),(tag,side,j)
   for key,value in expected.items():close(f[key],value,(tag,side,j,key))
   optional=['rawCellDispersion','cellDispersion','meanDonorMeasurementVarianceCPM','latentDonorRateVarianceCPM','newDonorRateVarianceCPM','gammaPriorShape','gammaPriorRatePerCPM'];assert all((k in f)==(k in expected) for k in optional)
   assert f['trainingCounts']==sum(g['counts'][j] for g in groups) and f['positiveDonors']==sum(g['counts'][j]>0 for g in groups)
  models[side]=model['features'];states[side]=dict(collections.Counter(f['status'] for f in model['features']))
 baseline=read(full/(fold['origin']+'-manifest.json.gz'))['genes'];changes=collections.Counter();phi_changes=[]
 for j,(c,t) in enumerate(zip(models['control'],models['treated'])):
  old=baseline[j];assert old['featureID']==c['featureID']==t['featureID'];a=old['controlCellDispersion'] is not None and old['treatedCellDispersion'] is not None;b='cellDispersion' in c and 'cellDispersion' in t;changes[f'full_{a}_fold_{b}']+=1
  for condition,f in [('control',c),('treated',t)]:
   x=old[condition+'CellDispersion'];y=f.get('cellDispersion')
   if x is not None and y is not None:phi_changes.append(abs(x-y)/max(1,abs(x)))
 folds.append(dict(fold,status='passed-independent-training-only-calibration',states=states,jointCellDispersionAvailability=dict(changes),comparableConditionGenes=len(phi_changes),dispersionChangesAbove1eMinus12=sum(x>1e-12 for x in phi_changes),maximumScaledDispersionChange=max(phi_changes)))
 print(tag,'verified',flush=True)
result={'status':'passed-all-training-only-calibrations','folds':folds,'geneConditionRecords':sum(2*f['featureCount'] for f in protocol['folds']),'numericalValuesCompared':values,'maximumScaledError':maximum,'seconds':time.time()-started,'jointLikelihoodFits':0,'biologicalOutcomeTests':0,'scope':'Training-only plug-in calibration inputs for donor-excluded joint fitting, not completed donor-held-out prediction or parameter uncertainty.'};(root/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='folds'}))
