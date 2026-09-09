#!/usr/bin/env python3
"""Diagnose the retained Kang NK gate failure without changing correction or gates.

Compare the original donor-held-out linear classifier with class-balanced linear
and exact 30-neighbor donor-held-out classifiers. All representations are frozen.
Source labels remain evaluation annotations, not authoritative inferred labels.
"""
import argparse,hashlib,json,platform
from importlib.metadata import version
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from check_full_integration import matrix
p=argparse.ArgumentParser(description=__doc__)
for name in ['reference','native','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
prior=json.loads((a.reference/'checks.json').read_text());inputs=np.load(a.reference/'evaluation-inputs.npz');x=inputs['scores'];donors=inputs['donors'];conditions=inputs['conditions'];types=inputs['cellTypes'];n,d=x.shape;levels=np.unique(types);truth=np.searchsorted(levels,types);nk=int(np.flatnonzero(levels=='NK cells')[0])
def save(name,v): (a.out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
spec=dict(schemaVersion=1,status='declared-before-diagnostic-results',priorReferenceSHA256=hashlib.sha256((a.reference/'checks.json').read_bytes()).hexdigest(),cells=n,components=d,seeds=[7,19,41],representations=['original']+[f'{kind}-{seed}' for kind in ['native','harmony'] for seed in [7,19,41]],
 classifiers=dict(original='StandardScaler + LogisticRegression C=1 lbfgs max_iter=2000; exactly reproduce original held-out donor predictions/recall.',balanced='Same fixed C, solver and scaling, with class_weight=balanced. Diagnostic sensitivity only; never replaces the original gate.',neighbors='Exact Euclidean 30-neighbor uniform vote over training donors, no scaling of latent axes, stable distance/source-index ties and lexicographic label ties.'),
 strata='Retain every donor, treatment and source cell type; no subset or threshold selected using diagnostic outputs.',localGeometry='Use every cell and the original retained global exact 30-neighbor lists. Summarize source-type composition in NK neighborhoods, by donor and condition.',
 acceptance='No new acceptance gate. Reproduce the original classifier metrics exactly; retain original failed gates even if a different diagnostic improves.',qualification='Exploratory failure localization on an already inspected cohort; not independent biological validation or an integration method change.')
save('design.json',spec)
representations={'original':x}
for seed in [7,19,41]:representations[f'native-{seed}']=matrix(a.native/f'seed-{seed}/scores.bin',n,d)
for seed in [7,19,41]:representations[f'harmony-{seed}']=np.load(a.reference/f'harmony-{seed}.npz')['scores']
summary=dict(status='in-progress',design=spec,versions={k:version(k) for k in ['numpy','scipy','scikit-learn']},platform=platform.platform(),types=levels.tolist(),runs=[])
for name,scores in representations.items():
 assert scores.shape==(n,d) and np.isfinite(scores).all();predictions={key:np.empty(n,dtype=np.int64) for key in ['original','balanced','neighbors']};fits=[]
 for donor in np.unique(donors):
  test=np.flatnonzero(donors==donor);train=np.flatnonzero(donors!=donor)
  assert len(np.unique(truth[train]))==len(levels)
  for kind,weight in [('original',None),('balanced','balanced')]:
   model=make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=2000,solver='lbfgs',class_weight=weight));model.fit(scores[train],truth[train]);assert model[-1].n_iter_.max()<2000
   predictions[kind][test]=model.predict(scores[test]);fits.append(dict(donor=str(donor),classifier=kind,trainingCells=len(train),testCells=len(test),iterations=int(model[-1].n_iter_.max()),means=model[0].mean_.tolist(),scales=model[0].scale_.tolist(),classes=model[-1].classes_.tolist(),coefficients=model[-1].coef_.tolist(),intercepts=model[-1].intercept_.tolist()))
  for start in range(0,len(test),64):
   selected=test[start:start+64];distance=cdist(scores[selected],scores[train])
   for offset,row in enumerate(distance):
    threshold=np.partition(row,29)[29];candidates=np.flatnonzero(row<=threshold);near=candidates[np.lexsort((train[candidates],row[candidates]))[:30]]
    predictions['neighbors'][selected[offset]]=np.bincount(truth[train[near]],minlength=len(levels)).argmax()
 np.savez_compressed(a.out/(name+'-predictions.npz'),truth=truth,donors=donors,conditions=conditions,**predictions);save(name+'-fits.json',fits)
 measurements={}
 for kind,prediction in predictions.items():
  folds=[];condition_folds=[]
  for donor in np.unique(donors):
   test=donors==donor;recalls={str(t):float(np.mean(prediction[test&(truth==i)]==i)) for i,t in enumerate(levels) if np.any(test&(truth==i))}
   confusion=np.bincount(truth[test]*len(levels)+prediction[test],minlength=len(levels)**2).reshape(len(levels),len(levels))
   folds.append(dict(donor=str(donor),cells=int(test.sum()),recalls=recalls,balancedAccuracy=float(np.mean(list(recalls.values()))),confusion=confusion.tolist()))
   for condition in np.unique(conditions):
    test=(donors==donor)&(conditions==condition)&(truth==nk)
    condition_folds.append(dict(donor=str(donor),condition=str(condition),cells=int(test.sum()),predictedTypeCounts={str(t):int(np.sum(prediction[test]==i)) for i,t in enumerate(levels)},nkRecall=float(np.mean(prediction[test]==nk)) if test.any() else None))
  measurements[kind]=dict(balancedAccuracy=float(np.mean([f['balancedAccuracy'] for f in folds])),recalls={str(t):float(np.mean([f['recalls'][str(t)] for f in folds if str(t) in f['recalls']])) for t in levels},folds=folds,nkConditionFolds=condition_folds)
 previous=prior['baseline'] if name=='original' else next(r['metrics'] for r in prior['native' if name.startswith('native') else 'references'] if r['seed']==int(name.rsplit('-',1)[1]))
 for t in levels:np.testing.assert_allclose(measurements['original']['recalls'][str(t)],previous['cellTypeRecall'][str(t)],rtol=0,atol=1e-15)
 np.testing.assert_allclose(measurements['original']['balancedAccuracy'],previous['cellTypeBalancedAccuracy'],rtol=0,atol=1e-15)
 label='baseline' if name=='original' else name
 neighbors=np.load(a.reference/(label+'-neighbors.npz'))['globalNeighbors'];assert neighbors.shape==(n,30)
 geometry=[]
 for donor in np.unique(donors):
  for condition in np.unique(conditions):
   mask=(donors==donor)&(conditions==condition)&(truth==nk);near=truth[neighbors[mask]]
   geometry.append(dict(donor=str(donor),condition=str(condition),cells=int(mask.sum()),neighborTypeFractions={str(t):float(np.mean(near==i)) for i,t in enumerate(levels)},meanDisplacement=float(np.mean(np.linalg.norm(scores[mask]-x[mask],axis=1)))))
 nk_neighbor_fraction=float(np.mean([np.mean(truth[neighbors[(donors==donor)&(truth==nk)]]==nk) for donor in np.unique(donors)]))
 run=dict(representation=name,classifiers=measurements,originalGateMetricsExactlyReproduced=True,nkNeighborhoods=geometry,donorBalancedNKNeighborFraction=nk_neighbor_fraction,scoreSHA256=hashlib.sha256(np.ascontiguousarray(scores,dtype='<f8').tobytes()).hexdigest())
 summary['runs'].append(run);save('checks.json',summary);print(json.dumps(dict(representation=name,nkRecall={k:v['recalls']['NK cells'] for k,v in measurements.items()},nkNeighborFraction=nk_neighbor_fraction)),flush=True)
summary['status']='frozen-representation-diagnosis-complete';summary['qualification']=spec['qualification'];save('checks.json',summary)
