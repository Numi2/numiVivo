"""Retrospective full-cohort uncertainty diagnosis from original sampled cells.

The widening diagnostic uses observed treated counts and is NOT a predictor.
"""
from pathlib import Path
import gzip, hashlib, json, os, time
import anndata as ad
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

root=Path(__file__).resolve().parent
study=Path('/Users/home/numivivo-gse181897-20260911')
out=root/'results';out.mkdir(exist_ok=False)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def write(name,d):
 (out/name).write_text(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+'\n')

protocol=json.loads((root/'protocol.json').read_bytes())
for p,h in protocol['sourceSHA256'].items():assert sha(study/p)==h,p
freeze=json.loads((study/'handoff/input-freeze.json').read_bytes())
for p,h in freeze['files'].items():assert sha(study/'handoff'/p)==h,p
pred_freeze=json.loads((study/'prediction-execution/prediction-freeze.json').read_bytes())
assert sha(study/'prediction-execution/retention.json')==pred_freeze['files']['retention.json']
retention=json.loads((study/'prediction-execution/retention.json').read_bytes())
def native(bundle,name):
 item=retention['bundles'][bundle][name]
 path=study/'prediction-execution/objects'/(item['SHA256']+'.gz')
 assert sha(path)==retention['objects'][item['SHA256']]['compressedSHA256']
 raw=gzip.decompress(path.read_bytes());assert hashlib.sha256(raw).hexdigest()==item['SHA256']
 return json.loads(raw)

cohort=json.loads((study/'prediction-inputs/cohort.json').read_bytes())
donors=cohort['eligibleDonors'];assert len(donors)==62
panel=json.loads((study/'prediction-inputs/panel.json').read_bytes());assert len(panel)==11800
groups=json.loads((study/'handoff/reference-groups.json').read_bytes())
lookup={(g['expID'],g['conditionCode']):i for i,g in enumerate(groups)}
reference=sparse.load_npz(study/'handoff/reference-pseudobulk.npz')
source=ad.read_h5ad(study/'handoff/B-source-codes-RNA.h5ad',backed='r')
assert source.shape==(15272,20303)
order=np.array([source.var_names.get_loc(x.removeprefix('symbol|')) for x in panel])
B=protocol['bootstrap']['replicates'];assert B==256
arrays={name:np.zeros((len(donors),len(panel))) for name in
 ['control','treated','controlVariance','treatedVariance','controlBootstrapMean','treatedBootstrapMean',
  'controlVarianceFirst','controlVarianceSecond','treatedVarianceFirst','treatedVarianceSecond']}
zero=np.zeros((62,11800,2),dtype=bool)
state=dict(status='running',pid=os.getpid(),startedUnix=time.time(),protocolSHA256=sha(root/'protocol.json'),
 scriptSHA256=sha(__file__),originalGroupsChecked=0,independentBootstrapCountChecks=0,
 predictionModelChanged=False,retrospectiveTreatedCountsUsed=True)
write('pipeline.json',state)
strata=[]
for di,donor in enumerate(donors):
 for ci,(code,label) in enumerate([('C','control'),('B','treated')]):
  group=groups[lookup[donor,code]];indices=np.array(group['sourceCellIndices'])
  assert len(indices)==group['cells'] and len(indices)>=2
  assert np.all(source.obs.iloc[indices]['exp_id'].astype(str).values==donor)
  assert np.all(source.obs.iloc[indices]['cond'].astype(str).values==code)
  X=source.X[indices,:].astype(np.int64).tocsr()
  sums=np.asarray(X.sum(axis=0,dtype=np.int64)).ravel()
  expected=reference[lookup[donor,code]].toarray().ravel()
  np.testing.assert_array_equal(sums,expected)
  assert int(sums.sum())==group['libraryCounts']
  original=np.log1p(sums[order]/sums.sum()*1e6)
  arrays[label][di]=original;zero[di,:,ci]=sums[order]==0
  rng=np.random.default_rng(np.random.SeedSequence([protocol['bootstrap']['seed'],int(donor),ci]))
  samples=np.zeros((B,len(panel)))
  n=len(indices)
  for begin in range(0,B,16):
   weights=rng.multinomial(n,np.full(n,1/n),size=min(16,B-begin))
   assert np.all(weights.sum(axis=1)==n)
   aggregate=np.asarray((X.T@weights.T).T,dtype=np.int64)
   totals=aggregate.sum(axis=1,dtype=np.int64)
   assert np.all(totals>0) and np.all(aggregate>=0)
   if begin==0:
    for k in range(2):
     expanded=np.repeat(np.arange(n),weights[k])
     np.testing.assert_array_equal(np.asarray(X[expanded].sum(axis=0)).ravel(),aggregate[k])
     state['independentBootstrapCountChecks']+=1
   samples[begin:begin+len(weights)]=np.log1p(aggregate[:,order]/totals[:,None]*1e6)
  arrays[label+'BootstrapMean'][di]=samples.mean(axis=0)
  arrays[label+'Variance'][di]=samples.var(axis=0,ddof=1)
  arrays[label+'VarianceFirst'][di]=samples[:B//2].var(axis=0,ddof=1)
  arrays[label+'VarianceSecond'][di]=samples[B//2:].var(axis=0,ddof=1)
  assert np.all(arrays[label+'Variance'][di]>=0)
  assert np.all(arrays[label+'Variance'][di][zero[di,:,ci]]==0)
  strata.append(dict(donor=donor,condition=label,cells=n,libraryCounts=int(sums.sum()),
     observedZeroFeatures=int(zero[di,:,ci].sum()),meanBootstrapVariance=float(arrays[label+'Variance'][di].mean())))
  state['originalGroupsChecked']+=1
 write('pipeline.json',state)
 if (di+1)%8==0: print(json.dumps({'donorsCompleted':di+1,'groupsChecked':state['originalGroupsChecked']}),flush=True)
source.file.close()
assert state['originalGroupsChecked']==124 and state['independentBootstrapCountChecks']==248
np.savez_compressed(out/'moments.npz',**arrays,zero=zero,donors=np.array(donors),panel=np.array(panel))
observed=arrays['treated']-arrays['control']
cell_variance=arrays['controlVariance']+arrays['treatedVariance']
records=[];summaries=[];categories=[]
for origin in ['Kang','HIRISA']:
 model=native('model-'+origin,'model.json')
 assert model['featureIDs']==panel
 mean=np.array(model['meanResponse']);variance=np.array([np.nan if v is None else v for v in model['donorResponseVariances']])
 n=len(model['trainingDonors']);available=np.isfinite(variance)
 first=native('prediction-'+origin+'-00','report.json')['predictions'][0]['meanResponsePredictiveInterval']
 critical=first['studentCriticalValue']
 base_variance=variance*(1+1/n)
 base_half=critical*np.sqrt(base_variance)
 residual=observed-mean
 mse=float(np.mean(residual**2));bias=float(np.mean(residual.mean(axis=0)**2));variation=float(np.mean(np.var(residual,axis=0)))
 assert abs(mse-bias-variation)<1e-12*max(1,mse)
 for di,donor in enumerate(donors):
  control,truth=arrays['control'][di],arrays['treated'][di]
  original=(np.abs(residual[di])<=base_half)&available
  widened=critical*np.sqrt(base_variance+cell_variance[di])
  raw=(np.abs(residual[di])<=widened)&available
  treated=(truth>=np.maximum(0,control+mean-widened))&(truth<=np.maximum(0,control+mean+widened))&available
  half_coverages=[]
  for half in ['First','Second']:
   sample_v=arrays['controlVariance'+half][di]+arrays['treatedVariance'+half][di]
   h=critical*np.sqrt(base_variance+sample_v)
   half_coverages.append(float(np.mean(np.abs(residual[di,available])<=h[available])))
  original_treated=(truth>=np.maximum(0,control+mean-base_half))&(truth<=np.maximum(0,control+mean+base_half))&available
  c=groups[lookup[donor,'C']]['cells'];t=groups[lookup[donor,'B']]['cells']
  records.append(dict(origin=origin,donor=donor,controlCells=c,treatedCells=t,availableFeatures=int(available.sum()),
    originalRawCoverage=float(original[available].mean()),originalTreatedCoverage=float(original_treated[available].mean()),
    observedCellWideningRawCoverage=float(raw[available].mean()),observedCellWideningTreatedCoverage=float(treated[available].mean()),
    observedCellWideningRawMeanWidth=float(np.mean(2*widened[available])),
    bootstrapHalvesRawCoverage=half_coverages,rawResponseRMSE=float(np.sqrt(np.mean(residual[di]**2))),
    meanObservedCellVariance=float(cell_variance[di].mean())))
 for name,mask in [('bothZero',zero[:,:,0]&zero[:,:,1]),('controlZeroOnly',zero[:,:,0]&~zero[:,:,1]),
                   ('treatedZeroOnly',~zero[:,:,0]&zero[:,:,1]),('bothNonzero',~zero[:,:,0]&~zero[:,:,1])]:
  measured=mask&available[None,:]
  categories.append(dict(origin=origin,stratum=name,donorFeaturePairs=int(mask.sum()),
    availablePairs=int(measured.sum()),rawResponseMSE=float(np.mean(residual[mask]**2)),
    originalRawCoverage=float(np.mean((np.abs(residual)<=base_half)[measured])),
    observedCellWideningRawCoverage=float(np.mean((np.abs(residual)<=critical*np.sqrt(base_variance+cell_variance))[measured])),
    meanObservedCellVariance=float(np.mean(cell_variance[mask]))))
 subset=[x for x in records if x['origin']==origin]
 old=json.loads((study/'prediction-scores-complete/summary.json').read_bytes())
 old=next(x for x in old if x['origin']==origin)
 for space,key in [('unclippedResponse','originalRawCoverage'),('predictedTreated','originalTreatedCoverage')]:
  np.testing.assert_allclose(np.mean([x[key] for x in subset]),old['intervalSummary'][space]['coverageMean'],rtol=0,atol=1e-14)
 summaries.append(dict(origin=origin,rawResponseMSE=mse,squaredDonorAverageResidual=bias,betweenDonorResidualVariance=variation,
   donorAverageResidualShare=bias/mse,meanObservedCellVariance=float(cell_variance.mean()),
   observedCellVarianceMagnitudeRelativeToMSE=float(cell_variance.mean()/mse),
   originalRawCoverage=float(np.mean([x['originalRawCoverage'] for x in subset])),
   observedCellWideningRawCoverage=float(np.mean([x['observedCellWideningRawCoverage'] for x in subset])),
   observedCellWideningTreatedCoverage=float(np.mean([x['observedCellWideningTreatedCoverage'] for x in subset])),
   widenedRawMeanWidth=float(np.mean([x['observedCellWideningRawMeanWidth'] for x in subset])),
   largestBootstrapHalfCoverageDifference=max(abs(x['bootstrapHalvesRawCoverage'][0]-x['bootstrapHalvesRawCoverage'][1]) for x in subset),
   cellCountVsRawRMSESpearman=float(spearmanr([min(x['controlCells'],x['treatedCells']) for x in subset],[x['rawResponseRMSE'] for x in subset]).statistic)))
write('donors.json',records);write('source-strata.json',strata);write('zero-strata.json',categories);write('summary.json',summaries)
state.update(status='completed-retrospective-cell-sampling-diagnostic',finishedUnix=time.time(),
    momentsSHA256=sha(out/'moments.npz'),originalIntervalCoverageReproduced=True,
    noNewBiologicalReplicates=True,noDeployableIntervalOrCalibrationClaim=True)
write('pipeline.json',state)
write('manifest.json',{'files':{p.name:dict(bytes=p.stat().st_size,SHA256=sha(p)) for p in sorted(out.iterdir()) if p.is_file()}})
print(json.dumps(summaries,indent=2),flush=True)
