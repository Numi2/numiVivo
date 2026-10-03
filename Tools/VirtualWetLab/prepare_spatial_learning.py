#!/usr/bin/env python3
"""Control-only spatial covariates for the existing native response learner.

Held-out query tensors contain no response RNA. Observation rows are materialized
only by the separately sealed evaluation path. Full cohort files remain immutable.
"""
import argparse,json
from pathlib import Path
import numpy as np
import anndata as ad
from sklearn.decomposition import PCA
from safetensors.numpy import save_file
from wetlab import read,write,sha,require

SEED=314159

def logcounts(a,rows):
 x=a.X[rows].toarray().astype(np.float32)
 return np.log1p(x/np.maximum(1,x.sum(1,keepdims=True))*1e6)

def prepare(cohort,out):
 cohort=Path(cohort);out=Path(out);out.mkdir(exist_ok=False)
 receipt=read(cohort/'cohort.json');protocol=read(cohort/'protocol.json')
 write(out/'preregistration.json',{'protocol':protocol,'cohortSHA256':sha(cohort/'cohort.json'),
   'instrumentRevision':sha(__file__),'clarifications':['same existing MLX architecture across ablations; Gaussian marginal variance not calibrated confidence',
   'all cell, reference and neighborhood inputs use admitted reference RNA only',
   'source presence mask applies to likelihood and metrics; no missing-as-zero evaluation',
   'training targets are anchor-level endpoint distributions; overlapping neighbors do not create independent replicates',
   'PCA uses training reference cells only; shared target embedding and role descriptor',
   'primary grouped evaluation uses unique outcome cell IDs within each target/type/neighborhood/role',
   'reference coverage is reported, not outcome-selected; inadequately covered populations cannot qualify',
   'cell-type-matched baseline backs off to role/target mean if training type absent',
   'external comparator is pinned State perturbation-mean; species-agnostic refit on this mouse axis, no transferred human weights'],
   'variantDescriptor':'reference cell PCA16 + reference neighborhood PCA16 + reference type probabilities7 + role1; frozen representation reduced to16 on training references when admitted',
   'distributions':'diagonal Gaussian in log1p(CPM), marginal samples do not imply joint gene or cell trajectories',
   'unseenTarget':'Cfap410; separate native run masks target identity and removes its training rows',
   'primarySelection':'equal-target validation mean whole-gene RMSE; three fixed checkpoints per variant',
   'testUse':'test observation extraction forbidden until prediction-seal.json exists'})
 targets=[x for x in receipt['targets'] if x['eligible']];targetindex={t['guideID']:i for i,t in enumerate(targets)}
 train=ad.read_h5ad(cohort/'chip2.h5ad');refs=np.flatnonzero(train.obs.reference_admitted.values)
 reference=logcounts(train,refs);variance=reference.var(0);panel=np.argsort(-variance,kind='stable')[:256];pca=PCA(n_components=16,svd_solver='full').fit(reference[:,panel]);scale=np.maximum(np.sqrt(pca.explained_variance_),.1)
 np.savez(out/'reference-pca.npz',panel=panel,mean=pca.mean_,components=pca.components_,scale=scale)
 features=train.var_names.tolist();write(out/'features.json',features)
 metricpanel=panel[:64].tolist();allmeta={};allarrays={};referenceexports={}
 del train
 for chip in ('chip2','chip3','chip1'):
  a=ad.read_h5ad(cohort/(chip+'.h5ad'));labels=a.obs.guide_assignment.astype(str).to_numpy();refs=np.flatnonzero(a.obs.reference_admitted.values);rlookup={int(r):i for i,r in enumerate(refs)};rv=logcounts(a,refs)
  require(np.array_equal(a.var_names,features),'Feature axis mismatch')
  refrows=np.vectorize(rlookup.__getitem__)(a.obsm['reference_indices']);weights=a.obsm['reference_weights'];emb=pca.transform(rv[:,panel])/scale
  local=(emb[refrows]*weights[:,:,None]).sum(1);neigh=a.obsm['neighbor_indices'];nhood=local[neigh].mean(1)
  referenceexports[chip]={'indices':refs.tolist(),'sourceSHA256':sha(cohort/(chip+'.h5ad'))}
  save_file({'counts':a.X[refs].toarray().astype(np.float32),'localIndices':refrows.astype(np.int32),'localWeights':weights.astype(np.float32),'neighbors':neigh.astype(np.int32)},str(out/(chip+'-reference.safetensors')))
  anchors=np.flatnonzero(np.isin(labels,list(targetindex)));meta=[];contexts=[];descriptors=[];tids=[];observed=[];obsvar=[]
  for cell in anchors:
   for role in ('direct','neighbor'):
    rows=np.array([cell]) if role=='direct' else neigh[cell][labels[neigh[cell]]=='barcode-negative']
    if not len(rows):continue
    # Both roles share reference state; only explicit neighborhood descriptor
    # is removed by the ablation. No target/neighbor expression enters here.
    context=(rv[refrows[cell]]*weights[cell,:,None]).sum(0)
    descriptor=np.r_[local[cell],nhood[cell],a.obsm['projected_type_probabilities'][cell],float(role=='neighbor')].astype(np.float32)
    refids=a.obsm['reference_indices'][cell].tolist();nrefids=np.unique(a.obsm['reference_indices'][neigh[cell]]).tolist()
    require(all(a.obs.reference_admitted.iloc[x] for x in refids+nrefids),'Non-reference input')
    # Barcode-negative neighbors of an intervention cannot also be admitted
    # reference cells, including reference neighbors shared across guides.
    require(not set(rows.tolist())&set(refids+nrefids),'Response enters its own reference input')
    contexts.append(context);descriptors.append(descriptor);tids.append(targetindex[labels[cell]])
    meta.append({'id':chip+':'+str(int(a.obs.source_cell_id.iloc[cell]))+':'+role,'chip':chip,'anchor':int(cell),'target':targets[tids[-1]]['target'],'targetIndex':tids[-1],'role':role,'cellType':str(a.obs.projected_cell_type.iloc[cell]),'annotationConfidence':float(a.obs.annotation_confidence.iloc[cell]),'referenceCoverage':bool(a.obs.reference_coverage.iloc[cell]),'neighborhood':int(nhood[cell,0]>0),'outcomeRows':rows.tolist(),'referenceRows':refids,'neighborhoodReferenceRows':nrefids,'controls':len(set(refids)),'xy':a.obsm['spatial'][cell].tolist()})
    if chip!='chip1':
     y=logcounts(a,rows);observed.append(y.mean(0));obsvar.append(y.var(0))
  context=np.asarray(contexts,np.float32);descriptor=np.asarray(descriptors,np.float32);tid=np.asarray(tids,np.int32);known=np.ones((len(tid),1),np.float32)
  arrays={'context':context,'descriptor':descriptor,'target':tid,'known':known}
  if chip!='chip1':arrays.update(observed=np.asarray(observed,np.float32),observedVariance=np.asarray(obsvar,np.float32),mask=np.tile(a.var.measured_in_source.to_numpy(np.float32),(len(tid),1)),stratum=(tid*2+descriptor[:,-1]).astype(np.int32))
  allarrays[chip]=arrays;allmeta[chip]=meta;write(out/(chip+'-rows.json'),meta)
  print(chip,len(meta),'anchor/role distributions',flush=True)
  del a,rv
 prior=np.asarray([np.mean((allarrays['chip2']['observed']-allarrays['chip2']['context'])[allarrays['chip2']['target']==i],axis=0) for i in range(len(targets))],np.float32)
 # Freeze target means and cell-type means using training observations only.
 baseline={'prior':prior};typekeys=sorted(set((m['target'],m['role'],m['cellType']) for m in allmeta['chip2']))
 write(out/'baseline-types.json',[list(x) for x in typekeys]);baseline['typePrior']=np.asarray([np.mean((allarrays['chip2']['observed']-allarrays['chip2']['context'])[[i for i,m in enumerate(allmeta['chip2']) if (m['target'],m['role'],m['cellType'])==key]],0) for key in typekeys],np.float32);save_file(baseline,str(out/'baselines.safetensors'))
 for variant in ('no-neighborhood','neighborhood','shuffled-neighborhood','unseen-Cfap410'):
  folder=out/variant;folder.mkdir();plan={'featureCount':len(features),'targetCount':len(targets),'descriptorCount':40,'hiddenWidth':64,'seed':SEED,'steps':[240,720,1440],'batchSize':16,'learningRate':.001,'weightDecay':.0001}
  write(folder/'plan.json',plan)
  for chip,raw in allarrays.items():
   arrays={k:v.copy() for k,v in raw.items()};d=arrays['descriptor']
   if variant=='no-neighborhood':d[:,16:32]=0
   if variant=='shuffled-neighborhood':
    rng=np.random.default_rng(SEED+int(chip[-1]));perm=np.arange(len(d))
    for role in (0,1):
     inds=np.flatnonzero(d[:,-1]==role);perm[inds]=rng.permutation(inds)
    d[:,16:32]=d[perm,16:32]
   if variant=='unseen-Cfap410':
    index=next(i for i,t in enumerate(targets) if t['target']=='Cfap410')
    if chip=='chip2':arrays={k:v[arrays['target']!=index] for k,v in arrays.items()}
    else:arrays['known'][arrays['target']==index]=0
   if chip=='chip2':
    arrays['prior']=prior.copy()
    if variant=='unseen-Cfap410':arrays['prior'][index]=0
   save_file(arrays,str(folder/(chip+'.safetensors')))
 write(out/'reference-rows.json',referenceexports)
 write(out/'prepared.json',{'format':'numivivo-spatial-learning-input/v1','cohort':str(cohort.resolve()),'cohortSHA256':sha(cohort/'cohort.json'),'protocolSHA256':sha(out/'preregistration.json'),'featuresSHA256':sha(out/'features.json'),'targets':targets,'metricPanelIndices':metricpanel,'rows':{k:len(v) for k,v in allmeta.items()},'files':{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cohort',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();prepare(a.cohort,a.output)
