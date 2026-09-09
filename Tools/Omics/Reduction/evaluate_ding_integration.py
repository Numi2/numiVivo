#!/usr/bin/env python3
"""Full-input, partial-source-label Ding evaluation with experiment-held-out readouts."""
import argparse,hashlib,json,platform,time
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path
import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from check_full_integration import matrix,independent

def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def exact_neighbors(x,experiments,types,known,k,workers):
 n=len(x);ids=np.arange(n);global_near=np.empty((n,k),dtype=np.int64);within=np.full((n,k),-1,dtype=np.int64)
 keys=np.array([json.dumps([e,t]) for e,t in zip(experiments,types)])
 groups={key:ids[known&(keys==key)] for key in np.unique(keys[known])}
 unavailable=[dict(stratum=json.loads(key),cells=len(rows),reason='Fewer than k+1 source-assigned cells; fixed k is not reduced') for key,rows in groups.items() if len(rows)<=k]
 def choose(row,candidates):
  values=row[candidates];threshold=np.partition(values,k-1)[k-1];tied=candidates[values<=threshold]
  return tied[np.lexsort((tied,row[tied]))[:k]]
 def block(start):
  dist=cdist(x[start:start+64],x)
  for offset,row in enumerate(dist):
   i=start+offset;row[i]=np.inf;global_near[i]=choose(row,ids)
   if known[i] and len(groups[keys[i]])>k:within[i]=choose(row,groups[keys[i]])
 with ThreadPoolExecutor(max_workers=workers) as pool:list(pool.map(block,range(0,n,64)))
 return global_near,within,groups,keys,unavailable

def metrics(x,inputs,k,workers,out,label):
 start=time.perf_counter();experiments=inputs['experiments'];methods=inputs['methods'];types=inputs['types'];known=inputs['known'];programs=inputs['programs'];program_names=inputs['programNames'];n=len(x)
 global_near,within,groups,keys,missing=exact_neighbors(x,experiments,types,known,k,workers)
 np.savez_compressed(out/(label+'-neighbors.npz'),globalNeighbors=global_near,stratifiedNeighbors=within)
 valid=within[:,0]>=0;excess=np.full(n,np.nan);entropy=np.full(n,np.nan);method_levels=np.unique(methods)
 for key,rows in groups.items():
  if len(rows)<=k:continue
  for i in rows:
   expected=(np.sum(methods[rows]==methods[i])-1)/(len(rows)-1);p=np.array([np.mean(methods[within[i]]==m) for m in method_levels])
   excess[i]=np.mean(methods[within[i]]==methods[i])-expected;entropy[i]=-np.sum(p[p>0]*np.log(p[p>0]))/np.log(len(method_levels))
 mixing_folds=[]
 for e in np.unique(experiments):
  for m in method_levels:
   select=valid&(experiments==e)&(methods==m)
   if select.any():mixing_folds.append(dict(experiment=str(e),method=str(m),cells=int(select.sum()),sameMethodExcess=float(np.mean(excess[select])),entropy=float(np.mean(entropy[select]))))
 def mean_mixing(key):return float(np.mean([np.mean([f[key] for f in mixing_folds if f['experiment']==e]) for e in np.unique(experiments)]))
 levels=np.unique(types[known]);truth=np.full(n,-1,dtype=np.int64);truth[known]=np.searchsorted(levels,types[known]);predictions={kind:np.full(n,-1,dtype=np.int64) for kind in ['original','balanced','neighbors']};fits=[];missing_classes=[]
 for experiment in np.unique(experiments):
  train=np.flatnonzero(known&(experiments!=experiment));test=np.flatnonzero(known&(experiments==experiment));assert len(train)>=k and len(test)>0
  absent=set(truth[test])-set(truth[train])
  missing_classes.extend(dict(experiment=str(experiment),cellType=str(levels[c]),reason='Source class absent in training experiments') for c in sorted(absent))
  for kind,weight in [('original',None),('balanced','balanced')]:
   model=make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=2000,solver='lbfgs',class_weight=weight));model.fit(x[train],truth[train]);assert model[-1].n_iter_.max()<2000
   predictions[kind][test]=model.predict(x[test]);fits.append(dict(experiment=str(experiment),classifier=kind,trainingCells=len(train),testCells=len(test),iterations=int(model[-1].n_iter_.max()),means=model[0].mean_.tolist(),scales=model[0].scale_.tolist(),classes=model[-1].classes_.tolist(),coefficients=model[-1].coef_.tolist(),intercepts=model[-1].intercept_.tolist()))
  def classify_block(start):
   selected=test[start:start+64];dist=cdist(x[selected],x[train])
   for offset,row in enumerate(dist):
    threshold=np.partition(row,k-1)[k-1];candidates=np.flatnonzero(row<=threshold);near=candidates[np.lexsort((train[candidates],row[candidates]))[:k]]
    predictions['neighbors'][selected[offset]]=np.bincount(truth[train[near]],minlength=len(levels)).argmax()
  with ThreadPoolExecutor(max_workers=workers) as pool:list(pool.map(classify_block,range(0,len(test),64)))
 classifiers={}
 for kind,prediction in predictions.items():
  folds=[]
  for experiment in np.unique(experiments):
   select=known&(experiments==experiment);recalls={str(t):float(np.mean(prediction[select&(truth==i)]==i)) for i,t in enumerate(levels) if np.any(select&(truth==i))}
   confusion=np.bincount(truth[select]*len(levels)+prediction[select],minlength=len(levels)**2).reshape(len(levels),len(levels))
   folds.append(dict(experiment=str(experiment),cells=int(select.sum()),balancedAccuracy=float(np.mean(list(recalls.values()))),recalls=recalls,confusion=confusion.tolist()))
  classifiers[kind]=dict(balancedAccuracy=float(np.mean([f['balancedAccuracy'] for f in folds])),recalls={str(t):float(np.mean([f['recalls'][str(t)] for f in folds if str(t) in f['recalls']])) for t in levels},folds=folds)
 np.savez_compressed(out/(label+'-predictions.npz'),truth=truth,experiments=experiments,known=known,**predictions);write(out/(label+'-fits.json'),fits)
 program_results={};missing_programs=[]
 def correlation(a,b,context):
  if len(a)<3 or len(np.unique(a))<2 or len(np.unique(b))<2:missing_programs.append(dict(**context,reason='Spearman needs >=3 cells and nonconstant measured and smoothed programs'));return None
  value=float(spearmanr(a,b).statistic);assert np.isfinite(value);return value
 for j,name in enumerate(program_names):
  values=programs[:,j];smoothed=values[global_near].mean(axis=1);by_experiment=[];by_stratum=[]
  for e in np.unique(experiments):
   select=experiments==e;c=correlation(values[select],smoothed[select],dict(program=str(name),experiment=str(e),scope='all-source-cells'))
   if c is not None:by_experiment.append(dict(experiment=str(e),spearman=c,cells=int(select.sum())))
  for key,rows in groups.items():
   if len(rows)<=k:continue
   c=correlation(values[rows],values[within[rows]].mean(axis=1),dict(program=str(name),stratum=json.loads(key),scope='source-assigned-stratum'))
   if c is not None:by_stratum.append(dict(stratum=json.loads(key),spearman=c,cells=len(rows)))
  program_results[str(name)]=dict(spearman=float(np.mean([v['spearman'] for v in by_experiment])) if by_experiment else None,withinStratumSpearman=float(np.mean([v['spearman'] for v in by_stratum])) if by_stratum else None,experiments=by_experiment,strata=by_stratum)
 nk=known&(types=='Natural killer cell');neighbor_known=known[global_near[nk]];nk_fraction=float(np.mean(types[global_near[nk]]=='Natural killer cell')) if nk.any() else None
 return dict(sameMethodExcess=mean_mixing('sameMethodExcess'),methodEntropy=mean_mixing('entropy'),mixingFolds=mixing_folds,cellTypeBalancedAccuracy=classifiers['original']['balancedAccuracy'],cellTypeRecall=classifiers['original']['recalls'],classifiers=classifiers,programs=program_results,nkSameTypeFractionAmongAllNeighbors=nk_fraction,nkNeighborSourceLabelCoverage=float(np.mean(neighbor_known)) if nk.any() else None,missingMixingStrata=missing,missingTrainingClasses=missing_classes,missingProgramStrata=missing_programs,evaluationSeconds=time.perf_counter()-start,executionWorkers=workers)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['prepare','native','reference'],required=True)
 for name in ['prepared','native','inputs','out']:p.add_argument('--'+name,type=Path,required=name in ['prepared','native','out'])
 p.add_argument('--workers',type=int,default=4);a=p.parse_args();assert 1<=a.workers<=8;a.out.mkdir(parents=True,exist_ok=False)
 protocol=json.loads((a.prepared/'protocol.json').read_text());n=protocol['cells'];k=protocol['evaluationNeighbors']
 versions={name:version(name) for name in ['numpy','scipy','scikit-learn','anndata','harmonypy']};assert versions['harmonypy']=='2.0.0'
 if a.stage=='prepare':
  artifacts=json.loads((a.native/'artifacts.json').read_text());assert any(Path(v['localPath'])==a.native/'pca' and v['allTransferredBytesExact'] for v in artifacts)
  source=ad.read_h5ad(a.prepared/'ding.h5ad');metadata=json.loads((a.native/'pca/metadata.json').read_text());d=json.loads((a.native/'pca/model.json').read_text())['options']['components']
  assert source.shape==(n,protocol['features']) and source.var_names.tolist()==[f['id'] for f in metadata['features']]
  assert [(str(s),str(b)) for s,b in zip(source.obs.native_sample,source.obs_names)]==[(c['sampleID'],c['barcode']) for c in metadata['cells']]
  assert all(c.get('group') is None for c in metadata['cells'])
  sample_map={s['id']:s for s in metadata['samples']}
  experiments=np.asarray(source.obs.experiment.astype(str),dtype=str);methods=np.asarray(source.obs.method.astype(str),dtype=str);types=np.asarray(source.obs.sourceCellType.astype(str),dtype=str);known=(source.obs.annotationStatus.astype(str)=='source-assigned').to_numpy()
  assert [sample_map[c['sampleID']]['condition'] for c in metadata['cells']]==experiments.tolist()
  assert [sample_map[c['sampleID']]['batchID'] for c in metadata['cells']]==methods.tolist()
  assert all(sample_map[c['sampleID']].get('donorID') is None for c in metadata['cells'])
  x=matrix(a.native/'pca/scores.bin',n,d);totals=np.asarray(source.X.sum(axis=1)).ravel();symbols=source.var.symbol.tolist();programs=[]
  for genes in protocol['programs'].values():
   assert all(symbols.count(g)==1 for g in genes);counts=source.X[:,[symbols.index(g) for g in genes]].toarray();programs.append(np.log1p(counts*(10000/totals[:,None])).mean(axis=1))
  inputs=dict(scores=x,experiments=experiments,methods=methods,types=types,known=known,programs=np.stack(programs,axis=1),programNames=np.array(list(protocol['programs'])))
  np.savez_compressed(a.out/'inputs.npz',**inputs)
  with np.load(a.out/'inputs.npz',allow_pickle=False) as reloaded:
   for key,value in inputs.items():
    assert not reloaded[key].dtype.hasobject
    np.testing.assert_array_equal(reloaded[key],value)
  baseline=metrics(x,inputs,k,a.workers,a.out,'baseline');write(a.out/'baseline.json',baseline);print('baseline measured',flush=True)
  erased=x.copy()
  for t in np.unique(types[known]):
   select=known&(types==t);erased[select]-=erased[select].mean(axis=0)
  control=metrics(erased,inputs,k,a.workers,a.out,'cellTypeErasure');write(a.out/'cellTypeErasure.json',control)
  write(a.out/'checks.json',dict(status='prepared',cells=n,sourceAssignedCells=int(known.sum()),sourceSHA256=sha(a.prepared/'ding.h5ad'),pcaScoresSHA256=sha(a.native/'pca/scores.bin'),metadataSHA256=sha(a.native/'pca/metadata.json'),protocolSHA256=sha(a.prepared/'protocol.json'),inputsSHA256=sha(a.out/'inputs.npz'),sourceAnnotationsExcludedFromNativeGroups=True,donorIdentityNotInvented=True,cellTypeErasureDetected=control['cellTypeBalancedAccuracy']<=baseline['cellTypeBalancedAccuracy']-protocol['gateMargins']['minimumNegativeControlAccuracyLoss'],versions=versions,platform=platform.platform()))
  print('erasure control measured',flush=True);return
 assert a.inputs is not None;checks=json.loads((a.inputs/'checks.json').read_text());assert checks['status']=='prepared' and sha(a.inputs/'inputs.npz')==checks['inputsSHA256'] and sha(a.prepared/'protocol.json')==checks['protocolSHA256']
 inputs=np.load(a.inputs/'inputs.npz');x=inputs['scores'];methods=inputs['methods'];baseline=json.loads((a.inputs/'baseline.json').read_text());margin=protocol['gateMargins']
 def gates(m):
  valid=lambda key:all(m['programs'][p][key] is not None and v[key] is not None and m['programs'][p][key]>=v[key]-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items())
  return dict(mixingImproved=m['sameMethodExcess']<baseline['sameMethodExcess'],cellTypesPreserved=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'],everyTypeRecallPreserved=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items()),programsPreserved=valid('spearman'),withinStratumProgramsPreserved=valid('withinStratumSpearman'))
 result=dict(status='in-progress',stage=a.stage,versions=versions,protocolSHA256=checks['protocolSHA256'],inputsSHA256=checks['inputsSHA256'],sourceAssignedCells=checks['sourceAssignedCells'],totalCells=n,cellTypeErasureDetected=checks['cellTypeErasureDetected'],runs=[])
 for mode,extra in protocol['modes'].items():
  for seed in protocol['seeds']:
   label=mode+'-'+str(seed);math=None;seconds=None
   if a.stage=='native':
    root=a.native/label;plan=json.loads((root/'plan.json').read_text())['integration'];assert plan==dict(protocol['integration'],**extra,seed=seed)
    assert json.loads((root/'metadata.json').read_text())==json.loads((a.native/'pca/metadata.json').read_text())
    y,math=independent(root,x,methods)
   else:
    opts=dict(theta=2,lamb=1 if mode=='fixed' else None,sigma=.1,nclust=100,tau=0,block_size=.05,max_iter_harmony=10,max_iter_kmeans=4,epsilon_cluster=.001,epsilon_harmony=.01,alpha=.2,batch_prop_cutoff=1e-5,ncores=1)
    start=time.perf_counter();h=harmonypy.run_harmony(x,pd.DataFrame({'method':methods}),['method'],**opts,random_state=seed,verbose=False);seconds=time.perf_counter()-start;y=np.asarray(h.Z_corr);assert y.shape==x.shape and np.isfinite(y).all();np.savez_compressed(a.out/(label+'-scores.npz'),scores=y)
    math=dict(objectives=list(h.objective_harmony))
   m=metrics(y,inputs,k,a.workers,a.out,label);result['runs'].append(dict(mode=mode,seed=seed,metrics=m,gates=gates(m),numericalReconstruction=math,referenceIntegrationSeconds=seconds));write(a.out/'checks.json',result);print(a.stage+' '+label+' measured',flush=True)
 result['status']='all-declared-runs-measured';result['allMeasuredAcceptanceGatesPassed']=all(all(r['gates'].values()) for r in result['runs']);result['completeSourceAnnotations']=checks['sourceAssignedCells']==n
 result['qualification']='Full UMI input; source labels available for only part of that input. Experiment-held-out classifiers and count-program preservation are separate readouts. The study does not establish donor identities or a treatment contrast. No prospective, fully annotated or cross-host speed qualification.';write(a.out/'checks.json',result)
if __name__=='__main__':main()
