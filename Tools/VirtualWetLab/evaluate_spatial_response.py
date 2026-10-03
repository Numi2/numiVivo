#!/usr/bin/env python3
"""Reveal source observations ONLY after every declared candidate is sealed."""
import argparse,math
from pathlib import Path
import numpy as np,anndata as ad
from scipy.stats import wasserstein_distance,norm
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,timestamp
from prepare_spatial_learning import logcounts

VARIANTS=['no-neighborhood','neighborhood','shuffled-neighborhood','frozen-Nicheformer','unseen-Cfap410']
def metrics(mean,var,observed,control,panel,mask):
 y=observed.mean(0);res=(mean-y)[mask];delta=(y-control)[mask];pd=(mean-control)[mask];sel=np.abs(delta)>.1
 out={'rmse':float(np.sqrt(np.mean(res**2))),'mae':float(np.mean(abs(res))),'effectSizeBias':float(res.mean()),'directionAgreement':float(np.mean(np.sign(pd[sel])==np.sign(delta[sel]))) if sel.any() else None,'directionGenes':int(sel.sum()),'genes':int(mask.sum())}
 if var is not None:
  v=np.maximum(var,1e-4);out['gaussianNLL']=float(.5*np.mean((np.log(2*np.pi*v)+(observed.var(0)+(mean-y)**2)/v)[mask]));out['marginal95Coverage']=float(np.mean((abs(observed[:,mask]-mean[mask])<=1.96*np.sqrt(v[mask]))));out['uncertaintyStatus']='uncalibrated marginal distribution; not animal-level confidence'
  q=norm.ppf((np.arange(64)+.5)/64);out['wassersteinPanelMean']=float(np.mean([wasserstein_distance(observed[:,g],mean[g]+np.sqrt(v[g])*q) for g in panel if mask[g]]))
 else:
  out['gaussianNLL']=None;out['marginal95Coverage']=None;out['wassersteinPanelMean']=float(np.mean([wasserstein_distance(observed[:,g],[mean[g]]) for g in panel if mask[g]]));out['uncertaintyStatus']='point baseline; no predictive variance'
 return out

def evaluate(root,out):
 root=Path(root);out=Path(out);out.mkdir(exist_ok=False);prepared=read(root/'prepared.json');cohort=Path(prepared['cohort']);require(sha(cohort/'cohort.json')==prepared['cohortSHA256'],'Cohort changed')
 rows=read(root/'chip1-rows.json');pred={};seals={}
 for variant in VARIANTS:
  folder=root/variant;seal=read(folder/'prediction-seal.json');require(sha(folder/'test-prediction/prediction.safetensors')==seal['predictionSHA256'],'Prediction changed');require(sha(Path(seal['weights']))==seal['weightsSHA256'],'Weights changed');require(sha(folder/'chip1.safetensors')==seal['inputsSHA256'],'Query changed');pred[variant]=load_file(str(folder/'test-prediction/prediction.safetensors'));seals[variant]=sha(folder/'prediction-seal.json')
 external=root/'external-State-PerturbMean';s=read(external/'prediction-seal.json');require(sha(external/'prediction.safetensors')==s['predictionSHA256'],'External prediction changed');pred['State-PerturbMean']=load_file(str(external/'prediction.safetensors'));seals['State-PerturbMean']=sha(external/'prediction-seal.json')
 query=load_file(str(root/'neighborhood/chip1.safetensors'));baseline=load_file(str(root/'baselines.safetensors'));keys=[tuple(x) for x in read(root/'baseline-types.json')];keyindex={x:i for i,x in enumerate(keys)}
 pred['no-change']={'mean':query['context']};pred['training-mean']={'mean':query['context']+baseline['prior'][query['target']]};matched=[]
 for i,r in enumerate(rows):
  ix=keyindex.get((r['target'],r['role'],r['cellType']));matched.append(query['context'][i]+(baseline['typePrior'][ix] if ix is not None else baseline['prior'][r['targetIndex']]))
 pred['cell-type-matched']={'mean':np.array(matched,np.float32)}
 baselinepath=out/'sealed-baselines.safetensors';save_file({k:v['mean'] for k,v in pred.items() if k in ['no-change','training-mean','cell-type-matched']},str(baselinepath))
 write(out/'prediction-seal.json',{'createdAt':timestamp(),'candidateSeals':seals,'baselinesSHA256':sha(baselinepath),'cohortSHA256':prepared['cohortSHA256'],'testOutcomesRead':False})
 # This is the first observation extraction. Queries and all model selections
 # above were already frozen. Preserve failures; no refitting below this line.
 a=ad.read_h5ad(cohort/'chip1.h5ad');observedRows=sorted(set(x for r in rows for x in r['outcomeRows']));y=logcounts(a,observedRows);yl={j:i for i,j in enumerate(observedRows)}
 trainmask=load_file(str(root/'neighborhood/chip2.safetensors'))['mask'][0].astype(bool);mask=a.var.measured_in_source.to_numpy()&trainmask;panel=prepared['metricPanelIndices'];groups={}
 for i,r in enumerate(rows):groups.setdefault((r['target'],r['cellType'],r['neighborhood'],r['role']),[]).append(i)
 summaries=[];arrays={};geneerrors=[]
 for number,(key,idx) in enumerate(sorted(groups.items())):
  outcomes=sorted(set(x for i in idx for x in rows[i]['outcomeRows']));obs=y[[yl[j] for j in outcomes]];control=query['context'][idx].mean(0);by={};population='|'.join(map(str,key));gid='population-'+str(number)
  for name,p in pred.items():
   if name=='unseen-Cfap410' and key[0]!='Cfap410':continue
   means=p['mean'][idx];m=means.mean(0);v=(p['variance'][idx]+means**2).mean(0)-m**2 if 'variance' in p else None
   by[name]=metrics(m,v,obs,control,panel,mask)
   if name=='neighborhood':arrays[gid+'-predicted']=m;arrays[gid+'-variance']=np.maximum(v,1e-4);arrays[gid+'-residual']=m-obs.mean(0);geneerrors.append((m-obs.mean(0))**2)
  arrays[gid+'-control']=control;arrays[gid+'-observed']=obs.mean(0)
  referenceIDs=set(x for i in idx for x in rows[i]['referenceRows']);covered=sum(rows[i]['referenceCoverage'] for i in idx)
  reasons=[]
  if by['no-neighborhood']['rmse']<by['neighborhood']['rmse']:reasons.append('neighborhood-free model has lower whole-gene RMSE')
  if by['cell-type-matched']['rmse']<by['neighborhood']['rmse']:reasons.append('cell-type-matched baseline has lower whole-gene RMSE')
  if by['shuffled-neighborhood']['rmse']<=by['neighborhood']['rmse']:reasons.append('neighborhood shuffling does not degrade RMSE')
  if len(outcomes)<5:reasons.append('fewer than five unique measured outcome cells')
  if covered<len(idx):reasons.append('some anchor reference distances exceed support limit')
  reasons.append('test chip has four isolated mSafe controls, below preregistered five')
  summaries.append({'id':gid,'target':key[0],'cellType':key[1],'neighborhood':key[2],'role':key[3],'chip':'chip1','animals':['paper-mouse1'],'outcomeCellCount':len(outcomes),'anchorCount':len(idx),'contributingReferenceCells':len(referenceIDs),'isolatedSafeControls':4,'referenceCoveredAnchors':covered,'previouslyExposed':key[0]=='Clu','metrics':by,'diagnostics':reasons,'support':'INSUFFICIENT CONTROL COVERAGE','biologicalPromotion':False,'outcomeRows':outcomes,'anchorRows':[rows[i]['anchor'] for i in idx]})
 save_file(arrays,str(out/'population-readouts.safetensors'));geneRMSE=np.sqrt(np.mean(geneerrors,0));genes=read(root/'features.json');order=np.argsort(-geneRMSE[mask]);valid=np.flatnonzero(mask)
 aggregate=[]
 for method in pred:
  included=[s for s in summaries if method in s['metrics']];targetScores=[]
  for target in sorted(set(s['target'] for s in included)):
   ss=[s for s in included if s['target']==target];targetScores.append({'target':target,'rmse':float(np.mean([s['metrics'][method]['rmse'] for s in ss])),'roleScores':{role:float(np.mean([s['metrics'][method]['rmse'] for s in ss if s['role']==role])) for role in ('direct','neighbor')}})
  aggregate.append({'method':method,'equalTargetRMSE':float(np.mean([s['rmse'] for s in targetScores])),'byTarget':targetScores,'track':'unseen training target, validation target exposed; exploratory specimen transfer only' if method=='unseen-Cfap410' else 'held-out specimen / previously trained interventions'})
 result={'format':'numivivo-spatial-response-evaluation/v1','predictionSealSHA256':sha(out/'prediction-seal.json'),'revealedAt':timestamp(),'observationsSHA256':sha(cohort/'chip1.h5ad'),'units':'log1p(CPM over each source measured gene universe)','commonScoredGenes':int(mask.sum()),'unmodeledGenes':int((~mask).sum()),'metricPanel':[genes[i] for i in panel],'populations':summaries,'aggregate':aggregate,'largestGeneErrors':[{'gene':genes[valid[i]],'rmse':float(geneRMSE[valid[i]])} for i in order[:30]],'evidence':'experimental candidate; no biological promotion','biologicalPromotion':False,'uncertainty':'uncalibrated; one held-out animal, no independent-animal confidence interval','v02Regression':'retained unchanged; not used to select this campaign','ridgeAdmission':'UNAVAILABLE on this split: legacy paired-donor ridge requires at least two independent training units; pooled chip2 cannot be relabelled as independent animals','spatialBenefitEstablished':False}
 write(out/'comparison.json',result);print([(r['method'],r['equalTargetRMSE']) for r in aggregate],flush=True)
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();evaluate(a.inputs,a.output)
