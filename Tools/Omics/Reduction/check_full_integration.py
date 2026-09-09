#!/usr/bin/env python3
"""Full-cohort donor integration evaluation with source labels used only for evaluation.

Fixed engineering preservation margins, not experimental causal validation.
The independent dense ridge solve is separate from the frozen trajectory oracle.
"""
import argparse,hashlib,json,platform,time
from importlib.metadata import version
from pathlib import Path
import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score,recall_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

def matrix(path,n,d):
 records=np.fromfile(path,dtype=[('row','<u4'),('column','<u4'),('value','<f8')]).reshape(n,d)
 assert np.equal(records['row'],np.arange(n)[:,None]).all() and np.equal(records['column'],np.arange(d)).all()
 assert np.isfinite(records['value']).all();return records['value'].copy()

def nearest(x,strata,k):
 n=len(x);global_indices=np.empty((n,k),dtype=np.int64);within=np.empty_like(global_indices);ids=np.arange(n)
 eligible={s:ids[strata==s] for s in np.unique(strata)}
 assert all(len(v)>k for v in eligible.values()),'Insufficient stratum for predeclared neighbor count'
 def choose(row,candidates):
  values=row[candidates];threshold=np.partition(values,k-1)[k-1];tied=candidates[values<=threshold]
  return tied[np.lexsort((tied,row[tied]))[:k]]
 for start in range(0,n,64):
  distances=cdist(x[start:start+64],x)
  for off,row in enumerate(distances):
   i=start+off;row[i]=np.inf
   global_indices[i]=choose(row,ids)
   within[i]=choose(row,eligible[strata[i]])
 return global_indices,within

def metrics(x,donors,conditions,types,program,k,out):
 levels=np.unique(donors);type_levels=np.unique(types);strata=np.array([json.dumps([t,c]) for t,c in zip(types,conditions)])
 global_neighbors,within=nearest(x,strata,k)
 np.savez_compressed(out,globalNeighbors=global_neighbors,stratifiedNeighbors=within)
 excess=[];entropy=[]
 for i,neighbors in enumerate(within):
  eligible=strata==strata[i];expected=((eligible&(donors==donors[i])).sum()-1)/(eligible.sum()-1)
  excess.append(np.mean(donors[neighbors]==donors[i])-expected)
  p=np.array([np.mean(donors[neighbors]==d) for d in levels]);entropy.append(-np.sum(p[p>0]*np.log(p[p>0]))/np.log(len(levels)))
 def donor_mean(v):return float(np.mean([np.mean(np.asarray(v)[donors==d]) for d in levels]))
 def classifier(train,test,target):
  model=make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=2000,solver='lbfgs'))
  model.fit(x[train],target[train]);assert model[-1].n_iter_.max()<2000,'Classifier did not converge'
  return model.predict(x[test]),int(model[-1].n_iter_.max())
 type_folds=[];condition_folds=[];missing=[]
 for donor in levels:
  test=donors==donor;train=~test
  if len(type_levels)>1:
   prediction,iterations=classifier(train,test,types)
   type_folds.append(dict(donor=str(donor),cells=int(test.sum()),balancedAccuracy=float(balanced_accuracy_score(types[test],prediction)),iterations=iterations,
    recalls={str(t):float(np.mean(prediction[types[test]==t]==t)) for t in type_levels if np.any(types[test]==t)}))
  for t in type_levels:
   tst=test&(types==t);trn=train&(types==t)
   if len(np.unique(conditions[tst]))!=2 or len(np.unique(conditions[trn]))!=2:
    missing.append(dict(donor=str(donor),cellType=str(t),reason='Treatment classifier needs both conditions in train and held-out donor'));continue
   prediction,iterations=classifier(trn,tst,conditions)
   condition_folds.append(dict(donor=str(donor),cellType=str(t),cells=int(tst.sum()),balancedAccuracy=float(balanced_accuracy_score(conditions[tst],prediction)),iterations=iterations))
 global_smoothed=program[global_neighbors].mean(axis=1);within_smoothed=program[within].mean(axis=1)
 global_corr=[];stratum_corr=[];unavailable=[]
 for donor in levels:
  select=donors==donor;global_corr.append(dict(donor=str(donor),spearman=float(spearmanr(program[select],global_smoothed[select]).statistic)))
  for t in type_levels:
   for c in np.unique(conditions):
    select=(donors==donor)&(types==t)&(conditions==c)
    if select.sum()<3 or len(np.unique(program[select]))<2 or len(np.unique(within_smoothed[select]))<2:
     unavailable.append(dict(donor=str(donor),cellType=str(t),condition=str(c),cells=int(select.sum()),reason='Spearman requires at least three cells and nonconstant program and prediction'));continue
    stratum_corr.append(dict(donor=str(donor),cellType=str(t),condition=str(c),cells=int(select.sum()),spearman=float(spearmanr(program[select],within_smoothed[select]).statistic)))
 by_type={str(t):float(np.mean([f['recalls'][str(t)] for f in type_folds if str(t) in f['recalls']])) for t in type_levels} if type_folds else {}
 condition_by_type={str(t):float(np.mean([f['balancedAccuracy'] for f in condition_folds if f['cellType']==str(t)])) for t in type_levels}
 return dict(sameDonorExcess=donor_mean(excess),donorEntropy=donor_mean(entropy),
  cellTypeBalancedAccuracy=float(np.mean([f['balancedAccuracy'] for f in type_folds])) if type_folds else None,cellTypeRecall=by_type,cellTypeFolds=type_folds,
  conditionBalancedAccuracy=float(np.mean(list(condition_by_type.values()))),conditionAccuracyByType=condition_by_type,conditionFolds=condition_folds,missingClassifierStrata=missing,
  programSpearman=float(np.mean([f['spearman'] for f in global_corr])),withinStratumProgramSpearman=float(np.mean([f['spearman'] for f in stratum_corr])),
  programDonorCorrelations=global_corr,programStratumCorrelations=stratum_corr,unavailableProgramStrata=unavailable,
  sameDonorExcessByType={str(t):float(np.mean([np.mean(np.asarray(excess)[(donors==d)&(types==t)]) for d in levels if np.any((donors==d)&(types==t))])) for t in type_levels})

def independent(root,x,donors):
 report=json.loads((root/'report.json').read_text());options=json.loads((root/'plan.json').read_text())['integration'];n,d=x.shape
 y=matrix(root/'scores.bin',n,d);r=matrix(root/'memberships.bin',n,options['clusters']);assignment=matrix(root/'assignment-scores.bin',n,d)
 batch=np.array(report['cellLevels']);assert np.asarray(report['levels'])[batch].tolist()==donors.tolist()
 np.testing.assert_allclose(r.sum(axis=1),1,rtol=0,atol=1e-12);assert (r>=0).all()
 sizes=np.bincount(batch);observed=np.stack([r[batch==b].sum(axis=0) for b in range(len(sizes))],axis=1);expected=r.sum(axis=0)[:,None]*sizes[None,:]/n
 distances=2*(1-assignment@np.asarray(report['assignmentCenters']).T)
 objective=(np.sum(r*distances)+options['temperature']*np.sum(r[r>0]*np.log(r[r>0]))+options['temperature']*options['diversity']*np.sum(observed*np.log((observed+expected+1)/(2*expected+1))))*2000/n
 np.testing.assert_allclose(objective,report['objectives'][-1],rtol=1e-10,atol=1e-10)
 reconstructed=x.copy()
 for c in range(r.shape[1]):
  active=np.flatnonzero(observed[c]/sizes>1e-5)
  if len(active)<2:continue
  masses=observed[c,active];sums=np.stack([r[batch==b,c]@x[batch==b] for b in active]);normal=np.diag(np.r_[masses.sum(),masses+options['ridge']]);normal[0,1:]=masses;normal[1:,0]=masses
  beta=np.linalg.solve(normal,np.vstack([sums.sum(axis=0),sums]))
  for slot,b in enumerate(active):
   mask=batch==b;reconstructed[mask]-=r[mask,c,None]*beta[slot+1]
 np.testing.assert_allclose(reconstructed,y,rtol=1e-10,atol=1e-10)
 improvements=-np.diff(report['objectives'])/np.maximum(np.abs(report['objectives'][:-1]),1e-12)
 np.testing.assert_allclose(improvements,report['relativeImprovements'],rtol=1e-12,atol=1e-12)
 if report['stoppingReason']=='relative-objective-tolerance':assert 0<=improvements[-1]<options['relativeTolerance']
 elif report['stoppingReason']=='objective-increase':assert improvements[-1]<0
 else:assert report['stoppingReason']=='iteration-limit' and len(improvements)==options['maximumIterations']
 return y,dict(independentObjective=float(objective),maximumCorrectionError=float(np.max(np.abs(reconstructed-y))),stoppingReason=report['stoppingReason'],objectives=report['objectives'])

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ['cohort','source','native','protocol','out']:p.add_argument('--'+name,type=Path if name!='cohort' else str,required=True)
 a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);protocol=json.loads(a.protocol.read_text());assert version('harmonypy')=='2.0.0'
 def save(v): (a.out/'checks.json').write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
 metadata=json.loads((a.native/'pca/metadata.json').read_text());model=json.loads((a.native/'pca/model.json').read_text());n=len(metadata['cells']);assert n==protocol['cohorts'][a.cohort]
 samples={s['id']:s for s in metadata['samples']};donors=np.array([samples[c['sampleID']]['donorID'] for c in metadata['cells']]);conditions=np.array([samples[c['sampleID']]['condition'] for c in metadata['cells']]);types=np.array([c.get('group','unreported') for c in metadata['cells']]);assert all(v not in ['','unreported'] for v in donors)
 x=matrix(a.native/'pca/scores.bin',n,model['options']['components']);source=ad.read_h5ad(a.source);assert source.n_obs==n
 mapping=json.loads((a.native/'pca/plan.json').read_text())['mapping'];barcodes=source.obs[mapping['barcodeColumn']].astype(str).tolist() if 'barcodeColumn' in mapping else source.obs_names.tolist()
 assert [(str(s),str(b)) for s,b in zip(source.obs[mapping['sampleColumn']],barcodes)]==[(c['sampleID'],c['barcode']) for c in metadata['cells']]
 assert source.var_names.tolist()==[f['id'] for f in metadata['features']]
 if 'groupColumn' in mapping:assert source.obs[mapping['groupColumn']].astype(str).tolist()==types.tolist()
 if a.cohort=='kang':assert source.obs['replicate'].astype(str).tolist()==donors.tolist() and source.obs['label'].astype(str).tolist()==conditions.tolist()
 genes=protocol['program'][a.cohort];assert all(source.var_names.tolist().count(g)==1 for g in genes)
 counts=sparse.csr_matrix(source.X);total=np.asarray(counts.sum(axis=1)).ravel();assert (total>0).all();selected=counts[:,[source.var_names.get_loc(g) for g in genes]].astype(float).toarray();program=np.log1p(selected*(10000/total[:,None])).mean(axis=1)
 np.savez_compressed(a.out/'evaluation-inputs.npz',scores=x,donors=donors,conditions=conditions,cellTypes=types,program=program)
 result=dict(status='in-progress',cohort=a.cohort,cells=n,dimensions=x.shape[1],protocolSHA256=hashlib.sha256(a.protocol.read_bytes()).hexdigest(),versions={s:version(s) for s in ['harmonypy','numpy','scipy','pandas','scikit-learn','anndata']},platform=platform.platform(),rareTypes=[str(t) for t in np.unique(types) if np.mean(types==t)<0.01],native=[],references=[])
 def measure(label,scores):
  start=time.perf_counter();m=metrics(scores,donors,conditions,types,program,protocol['evaluationNeighbors'],a.out/(label+'-neighbors.npz'));m['evaluationSeconds']=time.perf_counter()-start;print(a.cohort+' '+label+' measured',flush=True);return m
 baseline=measure('baseline',x);result['baseline']=baseline;save(result)
 for key,labels in [('conditionErasure',conditions),('cellTypeErasure',types)]:
  if key=='cellTypeErasure' and len(np.unique(types))==1:continue
  erased=x.copy()
  for label in np.unique(labels):
   mask=labels==label;erased[mask]-=erased[mask].mean(axis=0)
  result[key]=measure(key,erased);save(result)
 margin=protocol['gateMargins']
 def gates(m):
  values=dict(mixingImproved=m['sameDonorExcess']<baseline['sameDonorExcess'],conditionPreserved=m['conditionBalancedAccuracy']>=baseline['conditionBalancedAccuracy']-margin['maximumConditionBalancedAccuracyLoss'],programPreserved=m['programSpearman']>=baseline['programSpearman']-margin['maximumProgramSpearmanLoss'],withinStratumProgramPreserved=m['withinStratumProgramSpearman']>=baseline['withinStratumProgramSpearman']-margin['maximumProgramSpearmanLoss'],completeClassifierStrata=not m['missingClassifierStrata'])
  if baseline['cellTypeBalancedAccuracy'] is not None:
   values['cellTypesPreserved']=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'];values['everyTypeRecallPreserved']=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items())
  return values
 for seed in protocol['seeds']:
  y,math=independent(a.native/f'seed-{seed}',x,donors);m=measure(f'native-{seed}',y);result['native'].append(dict(seed=seed,metrics=m,gates=gates(m),independentReconstruction=math));save(result)
 opts=dict(theta=2,lamb=1,sigma=0.1,nclust=100,tau=0,block_size=0.05,max_iter_harmony=10,max_iter_kmeans=4,epsilon_cluster=0.001,epsilon_harmony=0.01,alpha=0.2,batch_prop_cutoff=1e-5,ncores=1)
 for seed in protocol['seeds']:
  start=time.perf_counter();h=harmonypy.run_harmony(x,pd.DataFrame({'donor':donors}),['donor'],**opts,random_state=seed,verbose=False);seconds=time.perf_counter()-start;y=np.asarray(h.Z_corr);assert y.shape==x.shape and np.isfinite(y).all();np.savez_compressed(a.out/f'harmony-{seed}.npz',scores=y)
  m=measure(f'harmony-{seed}',y);result['references'].append(dict(seed=seed,metrics=m,gates=gates(m),integrationSeconds=seconds,objectives=list(h.objective_harmony)));save(result)
 result['negativeControlsDetectLoss']=dict(condition=result['conditionErasure']['conditionBalancedAccuracy']<baseline['conditionBalancedAccuracy']-margin['minimumNegativeControlAccuracyLoss'])
 if 'cellTypeErasure' in result:result['negativeControlsDetectLoss']['cellType']=result['cellTypeErasure']['cellTypeBalancedAccuracy']<baseline['cellTypeBalancedAccuracy']-margin['minimumNegativeControlAccuracyLoss']
 result['status']='full-cohort-native-and-reference-measured';result['allNativePreservationGatesPassed']=all(all(r['gates'].values()) for r in result['native'])
 result['qualification']='Three native and three independent Harmony initializations. Native numerical reconstruction and biological preservation are separate gates. Transductive correction sees all unlabeled cells; held-out donor classifiers do not establish prospective prediction. Source cell labels are evaluation annotations, not newly authoritative inferred labels. Missing/constant program strata remain unavailable. CPU reference and native benchmark ran on different hosts; no speed comparison.';save(result)
 print(json.dumps(dict(cohort=a.cohort,allNativePreservationGatesPassed=result['allNativePreservationGatesPassed'],negativeControls=result['negativeControlsDetectLoss'])),flush=True)
if __name__=='__main__':main()
