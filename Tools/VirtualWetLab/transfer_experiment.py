#!/usr/bin/env python3
"""One nested, grouped technical-transfer experiment over exposed development data.

Frozen ESM2 versus coarse descriptors and target IDs, same native owner and
matched ridge information. No new architecture, optimizer search or independent
outcomes. All selectors precede outer scoring; target/guide/bag groups never split.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,timestamp
from prepare_intervention_design import load_source,axis,normalized
from train_intervention_design import native
from train_cohort_recovery import measure,predict
from design_objective import objective,utility,selection_score
SEED=271828
STEPS=[240,1440,5760]
REPRESENTATIONS=['coarse','ESM2','target-ID']
def key(g):return (g['source'],g['context'],g['target'])
def context(g):return g['source']+'|'+g['context']
def unique_rows(rows):return [next(r for r in rows if key(r)==k) for k in sorted(set(map(key,rows)))]
def stable(xs):return sorted(xs,key=lambda s:hashlib.sha256(('transfer-v1:'+s).encode()).hexdigest())

def register(inputs,embeddings,targets,binary,out):
 require(not out.exists(),'Retain previous campaign');out.mkdir()
 groups=unique_rows(read(inputs/'training-rows.json'));targetsList=sorted({g['target'] for g in groups});folds=[]
 order=stable(targetsList)
 for f in range(3):
  test=[i for i,g in enumerate(groups) if order.index(g['target'])%3==f];train=[i for i in range(len(groups)) if i not in test];valTargets=set(stable({groups[i]['target'] for i in train})[::4]);val=[i for i in train if groups[i]['target'] in valTargets];fit=[i for i in train if i not in val]
  folds.append({'id':'held-target-'+str(f),'task':'unseen-target-familiar-context','training':train,'innerTraining':fit,'innerValidation':val,'test':test,'grouping':'target globally, all contexts, guides, bags and source rows together'})
 exclusions=[]
 for ctx in sorted(set(map(context,groups))):
  testAll=[i for i,g in enumerate(groups) if context(g)==ctx];train=[i for i,g in enumerate(groups) if context(g)!=ctx];known={groups[i]['target'] for i in train};test=[i for i in testAll if groups[i]['target'] in known]
  if len(test)<5:exclusions.append({'context':ctx,'knownTargetGroups':len(test),'reason':'Fewer than five overlapping target groups; no adequate known-target context test'});continue
  innerContext='GSE90063|K562-high-MOI';fit=[i for i in train if context(groups[i])!=innerContext];knownInner={groups[i]['target'] for i in fit};val=[i for i in train if context(groups[i])==innerContext and groups[i]['target'] in knownInner]
  require(len(val)>0,'No nested context selector');folds.append({'id':'held-context-'+ctx.split('|')[-1],'task':'known-target-unseen-context','training':train,'innerTraining':fit,'innerValidation':val,'test':test,'excludedTestGroups':sorted(set(testAll)-set(test)),'grouping':'complete source/context, all guides and bags','selectorLimit':'Only three shared targets in inner K562-high-MOI context; mixed study/modality shift. Not broadly representative context validation.'})
 assert len(groups)==81
 for f in folds:
  assert not set(f['training'])&set(f['test']) and not set(f['innerTraining'])&set(f['innerValidation'])
  if f['task'].startswith('unseen-target'):assert not {groups[i]['target'] for i in f['training']}&{groups[i]['target'] for i in f['test']}
  else:assert not {context(groups[i]) for i in f['training']}&{context(groups[i]) for i in f['test']}
 specs=[s for s in read(inputs/'prepared.json')['sources'] if s['split']=='training']
 write(out/'registration.json',{'createdAt':timestamp(),'question':'Does frozen protein sequence information transfer to unseen targets separately from known-target context transfer?','scope':'Nested technical development transfer on previously exposed data; not independent biological validation','hypothesis':'ESM2 lowers held-target response RMSE versus coarse and target-ID controls, retains direction, and loses performance when target inputs are exchanged. Context transfer assessed separately.','refutation':'No lower error than no change and matched coarse/ridge controls, no direction improvement, or correct target no better than wrong target: no transferable benefit for that task. Selection gains remain separate outcomes.','binary':str(binary),'binarySHA256':sha(binary),'inputs':str(inputs),'inputRegistrationSHA256':sha(inputs/'prepared.json'),'targetManifest':str(targets),'targetManifestSHA256':sha(targets),'embeddings':str(embeddings),'embeddingsSHA256':sha(embeddings/'embeddings.safetensors'),'embeddingProvenanceSHA256':sha(embeddings/'provenance.json'),'ownerSHA256':sha(__file__),'groups':groups,'sources':specs,'folds':folds,'contextExclusions':exclusions,'representations':REPRESENTATIONS,'linearAlphas':[.1,1.,10.],'native':{'hiddenWidth':64,'learningRate':.001,'weightDecay':.0001,'optimizer':'adam','objective':'mean','sampling':'group','batchSize':16,'steps':STEPS},'preprocessing':'Common gene identities from source metadata; 512 highest-variance training-context control genes plus original immediate-early objective. Fit feature selection, response RMS scale and ESM PCA16 separately within inner-training and refit outer-training; unique training targets only for PCA. No held-out perturbed RNA in input context.','objective':{'genes':['FOS','JUN','JUNB','EGR1'],'preserveGenes':[],'penalty':0},'selection':'Separate per fold and representation; inner equal-group response RMSE, ties earlier step or stronger ridge penalty. Refit selected settings on outer training before sealing predictions. Never select using outer outcomes or combine task scores.','metrics':'Original group response RMSE, endpoint RMSE, direction deadband0.1, nearest response target within source/context, correct versus exchanged-target input, original selection_score utility/regret/random/no-intervention/training-mean ranking. No response-variation selection.','uncertainty':'No biological confidence intervals; resampled bags/cells/genes are not biological units. IDs remain source-bound with unresolved replicate independence.','budget':'5 folds x 3 native representations: one inner fit through 5760 and one outer refit at selected checkpoint; 3 ridge penalties per representation inside each fold. No architecture/optimizer sweep or restarts.','promotion':'No biological promotion; freeze selectors and report both tasks before reserving any new independent evaluation'})

def prepare(out):
 reg=read(out/'registration.json');require(sha(__file__)==reg['ownerSHA256'],'Registered owner changed');dest=out/'raw';require(not dest.exists(),'Retain raw preparation');dest.mkdir();groups=reg['groups'];axes=[]
 for s in reg['sources']:
  require(sha(s['path'])==s['sha256'],'Source changed');a=load_source(s);axes.append(set(axis(a,{},s['species'])));del a
 genes=sorted(set.intersection(*axes));write(dest/'features.json',genes)
 bagY=[];bagVar=[];base=[];controlVar={};groupMeta=[];observed=[];coverage=[];bagProvenance=[]
 for s in reg['sources']:
  a=load_source(s);ax=axis(a,{},s['species']);cols=[ax[g] for g in genes]
  for gi,g in [(i,g) for i,g in enumerate(groups) if g['source']==s['id']]:
   ckey=context(g)
   if ckey not in controlVar:
    c=normalized(a,g['controls'],cols);controlVar[ckey]=c.var(0);controlMean=c.mean(0);detected=(c>0).sum(0);del c
   rng=np.random.default_rng(SEED+gi);indices=[rng.choice(g['rows'],min(32,len(g['rows'])),replace=False) for _ in range(8)];ys=[normalized(a,x,cols) for x in indices];bagY.append(np.asarray([x.mean(0) for x in ys]));bagVar.append(np.asarray([x.var(0) for x in ys]));base.append(controlMean);coverage.append(detected)
   # Population endpoint is separate from eight overlapping training bags.
   total=np.zeros(len(genes),np.float64)
   for j in range(0,len(g['rows']),128):total+=normalized(a,g['rows'][j:j+128],cols).sum(0)
   observed.append(total/len(g['rows']));groupMeta.append(gi)
   guideColumn='perturbation' if 'perturbation' in a.obs else 'target';guides=sorted(set(a.obs.iloc[g['rows']][guideColumn].astype(str)));bagProvenance.append({'groupIndex':gi,'sourceSHA256':s['sha256'],'guides':guides,'sourceRows':g['rows'],'bagSourceRows':[x.tolist() for x in indices],'unit':g['unit'],'scope':'source row identities; technical sample, no independent biological unit claim'})
  del a;print(s['id'],'prepared',flush=True)
 order=np.argsort(groupMeta);save_file({'bags':np.asarray(bagY,np.float32)[order],'bagVariance':np.asarray(bagVar,np.float32)[order],'context':np.asarray(base,np.float32)[order],'observed':np.asarray(observed,np.float32)[order],'detectedControls':np.asarray(coverage,np.int32)[order]},str(dest/'population.safetensors'));save_file(controlVar,str(dest/'control-variance.safetensors'));write(dest/'source-groups.json',sorted(bagProvenance,key=lambda x:x['groupIndex']));write(dest/'seal.json',{'files':{p.name:sha(p) for p in dest.iterdir() if p.is_file()},'scope':'Exposed development population observations; outer folds excluded from all fitting'})

def features_for(groups,ids,features,variances,objectiveGenes):
 contexts=sorted({context(groups[i]) for i in ids});score=sum(variances[x] for x in contexts);order=sorted(range(len(features)),key=lambda i:(-float(score[i]),features[i]));chosen=order[:512]
 for gene in objectiveGenes:
  require(gene in features,'Declared objective unavailable; never substitute');j=features.index(gene)
  if j not in chosen:chosen.append(j)
 return np.asarray(chosen,int)

def preprocessing(reg,raw,features,var,ids,rep,embeddings,folder):
 groups=reg['groups'];targetNames=sorted(read(reg['targetManifest']));sourceTargets=read(reg['targetManifest']);trainTargets=sorted({groups[i]['target'] for i in ids});cols=features_for(groups,ids,features,var,reg['objective']['genes']);bio={t:np.zeros(149,np.float32) for t in targetNames};provenance={'trainingGroups':ids,'trainingTargets':trainTargets,'features':[features[j] for j in cols],'representation':rep}
 if rep=='coarse':bio={t:np.asarray(sourceTargets[t]['descriptor'],np.float32) for t in targetNames}
 elif rep=='ESM2':
  x=np.asarray([embeddings[t] for t in trainTargets]);center=x.mean(0);u,s,v=np.linalg.svd(x-center,full_matrices=False);n=min(16,len(trainTargets)-1);components=v[:n];coords=(x-center)@components.T;scale=max(float(np.sqrt(np.mean(coords**2))),1e-6)
  for t in targetNames:bio[t][:n]=((embeddings[t]-center)@components.T)/scale
  save_file({'center':center,'components':components,'scale':np.asarray([scale],np.float32)},str(folder/'pca.safetensors'));provenance.update(pcaTrainingTargets=trainTargets,pcaComponents=n,pcaSHA256=sha(folder/'pca.safetensors'))
 y=raw['bags'][ids][:,:,cols]-raw['context'][ids][:,None,cols];responseScale=np.maximum(.1,np.sqrt(np.mean(y*y,axis=(0,1)))).astype(np.float32);save_file({'responseScale':responseScale},str(folder/'scaling.safetensors'));write(folder/'preprocessing.json',provenance)
 return cols,bio,targetNames,trainTargets,responseScale

def data_for(reg,raw,ids,prep,rep,bags=False,wrong=False):
 cols,bio,targets,trainTargets,scale=prep;groups=reg['groups'];records=[];a={k:[] for k in ['context','descriptor','target','known','observed','observedVariance','mask','stratum']};swap={i:i for i in ids}
 if wrong:
  for ctx in sorted({context(groups[i]) for i in ids}):
   ii=sorted([i for i in ids if context(groups[i])==ctx],key=lambda i:groups[i]['target'])
   if len(ii)>1:swap.update(dict(zip(ii,ii[1:]+ii[:1])))
 for local,i in enumerate(ids):
  g=groups[i];target=groups[swap[i]]['target'];meta=[float(g['modality']=='CRISPR-KO'),float(g['modality']=='CRISPRi'),0.,float(g['species']=='mouse'),float(g['species']=='human'),float(g['context']=='stimulated'),float(g['context']=='unstimulated')]
  for b in range(8 if bags else 1):
   a['context'].append(raw['context'][i,cols]);a['descriptor'].append(np.r_[bio[target],meta]);a['target'].append(targets.index(target));a['known'].append([float(rep=='target-ID' and target in trainTargets)]);a['stratum'].append(local);a['observed'].append(raw['bags'][i,b,cols] if bags else raw['observed'][i,cols]);a['observedVariance'].append(raw['bagVariance'][i,b,cols] if bags else np.zeros(len(cols)));a['mask'].append(np.ones(len(cols)));records.append({k:g[k] for k in ['source','context','target','unit']})
 a={k:np.asarray(v,np.int32 if k in ['target','stratum'] else np.float32) for k,v in a.items()};a['prior']=np.zeros((len(targets),len(cols)),np.float32);a['responseScale']=scale;return a,records

def query(a):return {k:a[k] for k in ['context','descriptor','target','known']}
def metric(pred,a,rows):
 result=measure(pred,a,rows);result['endpointRMSE']=float(np.sqrt(np.mean((pred-a['observed'])**2)));return result

def linearX(a,rep,count):
 block=np.eye(count,dtype=np.float32)[a['target']]*a['known'] if rep=='target-ID' else a['descriptor'][:,:149]
 return np.c_[a['context']/10,block,a['descriptor'][:,149:]].astype(np.float64)
def ridge_fit(x,y,alpha):
 xm=x.mean(0);ym=y.mean(0);x=x-xm;y=y-ym;w=x.T@np.linalg.solve(x@x.T+alpha*np.eye(len(x)),y);return xm,ym,w

def decision(reg,features,ids,mean,a,trainIds,raw,cols):
 groups=reg['groups'];idx=objective(features,reg['objective']);p=utility(mean-a['context'],None,idx)[0];y=utility(a['observed']-a['context'],None,idx)[0];scores=[]
 for ctx in sorted({context(groups[i]) for i in ids}):
  ii=[j for j,i in enumerate(ids) if context(groups[i])==ctx];candidates=[groups[ids[j]]['target'] for j in ii]+['no-intervention'];baseline=[]
  for j in ii:
   g=groups[ids[j]];match=[k for k in trainIds if groups[k]['target']==g['target'] and groups[k]['modality']==g['modality']];scope='same target and modality' if match else 'training modality mean'
   if not match:match=[k for k in trainIds if groups[k]['modality']==g['modality']]
   delta=(raw['observed'][match][:,cols]-raw['context'][match][:,cols]).mean(0);baseline.append(float(utility(delta,None,idx)[0]))
  score=selection_score(candidates,list(p[ii])+[0.],list(y[ii])+[0.],baseline+[0.]);score.update(context=ctx,candidates=candidates,controlCounts=[len(groups[ids[j]]['controls']) for j in ii],outcomeCounts=[len(groups[ids[j]]['rows']) for j in ii],trainingBaseline='same target and modality where admitted; otherwise training modality mean');scores.append(score)
 return scores

def execute(out):
 reg=read(out/'registration.json');require(sha(__file__)==reg['ownerSHA256'],'Owner changed');binary=Path(reg['binary']);require(sha(binary)==reg['binarySHA256'],'Native runtime changed');emb=load_file(str(Path(reg['embeddings'])/'embeddings.safetensors'));require(sha(Path(reg['embeddings'])/'embeddings.safetensors')==reg['embeddingsSHA256'],'Embedding changed');raw=load_file(str(out/'raw/population.safetensors'));features=read(out/'raw/features.json');var=load_file(str(out/'raw/control-variance.safetensors'));results=[]
 for fold in reg['folds']:
  fd=out/fold['id'];fd.mkdir(exist_ok=True)
  for rep in reg['representations']:
   dest=fd/rep
   if (dest/'result.json').exists():results.append(read(dest/'result.json'));continue
   require(not dest.exists(),'Interrupted model retained; explicit recovery required');dest.mkdir();inner=dest/'inner';inner.mkdir();prep=preprocessing(reg,raw,features,var,fold['innerTraining'],rep,emb,inner)
   train,trainRows=data_for(reg,raw,fold['innerTraining'],prep,rep,True);val,valRows=data_for(reg,raw,fold['innerValidation'],prep,rep);plan={**reg['native'],'featureCount':len(prep[0]),'targetCount':len(prep[2]),'descriptorCount':156,'seed':SEED,'trainingBudget':STEPS[-1],'diagnostics':False};write(inner/'plan.json',plan);save_file(train,str(inner/'training.safetensors'));save_file(query(val),str(inner/'validation-query.safetensors'));native(binary,'train',inner/'plan.json',inner/'training.safetensors',inner/'training');checks=[]
   for step in STEPS:
    pred,digest=predict(binary,inner/'plan.json',inner/'validation-query.safetensors',inner/'training'/f'weights-{step}.safetensors');checks.append({'step':step,'validation':metric(pred['mean'],val,valRows),'predictionSHA256':digest})
   selected=min(checks,key=lambda c:(c['validation']['equalGroupRMSE'],c['step']));small,_=data_for(reg,raw,fold['innerTraining'],prep,rep);x=linearX(small,rep,len(prep[2]));vx=linearX(val,rep,len(prep[2]));linear=[]
   for alpha in reg['linearAlphas']:
    center,intercept,w=ridge_fit(x,small['observed']-small['context'],alpha);pred=val['context']+(vx-center)@w+intercept;linear.append({'alpha':alpha,'validation':metric(pred.astype(np.float32),val,valRows)})
   chosen=min(linear,key=lambda c:(c['validation']['equalGroupRMSE'],-c['alpha']));write(dest/'selector.json',{'frozenAt':timestamp(),'outerOutcomesUsed':False,'native':selected,'nativeCheckpoints':checks,'linear':chosen,'linearCandidates':linear,'selectionScope':fold.get('selectorLimit',fold['task'])})
   final=dest/'refit';final.mkdir();prep=preprocessing(reg,raw,features,var,fold['training'],rep,emb,final);train,trainRows=data_for(reg,raw,fold['training'],prep,rep,True);test,testRows=data_for(reg,raw,fold['test'],prep,rep);wrong,_=data_for(reg,raw,fold['test'],prep,rep,wrong=True);plan.update(featureCount=len(prep[0]),steps=[selected['step']],trainingBudget=selected['step']);write(final/'plan.json',plan);save_file(train,str(final/'training.safetensors'));save_file(query(test),str(final/'test-query.safetensors'));save_file(query(wrong),str(final/'wrong-target-query.safetensors'));native(binary,'train',final/'plan.json',final/'training.safetensors',final/'training');weights=final/'training'/f"weights-{selected['step']}.safetensors";native(binary,'predict',final/'plan.json',final/'test-query.safetensors',final/'prediction',weights);native(binary,'predict',final/'plan.json',final/'wrong-target-query.safetensors',final/'wrong-target-prediction',weights);small,_=data_for(reg,raw,fold['training'],prep,rep);center,intercept,w=ridge_fit(linearX(small,rep,len(prep[2])),small['observed']-small['context'],chosen['alpha']);save_file({'center':center,'intercept':intercept,'weights':w},str(final/'linear.safetensors'))
   lp=test['context']+(linearX(test,rep,len(prep[2]))-center)@w+intercept;wp=test['context']+(linearX(wrong,rep,len(prep[2]))-center)@w+intercept;save_file({'mean':lp.astype(np.float32),'wrongTarget':wp.astype(np.float32)},str(final/'linear-predictions.safetensors'));write(dest/'prediction-seal.json',{'frozenAt':timestamp(),'selectorSHA256':sha(dest/'selector.json'),'files':{str(p.relative_to(final)):sha(p) for p in final.rglob('*') if p.is_file()},'scope':'development prediction sealed before outer evaluation'});
   models={'native':load_file(str(final/'prediction/prediction.safetensors'))['mean'],'ridge':lp.astype(np.float32)};wrongs={'native':load_file(str(final/'wrong-target-prediction/prediction.safetensors'))['mean'],'ridge':wp.astype(np.float32)};result={'fold':fold['id'],'task':fold['task'],'representation':rep,'selectedStep':selected['step'],'selectedAlpha':chosen['alpha'],'featureCount':len(prep[0]),'testGroups':len(fold['test']),'models':{},'biologicalPromotion':False}
   for name,p in models.items():result['models'][name]={'metrics':metric(p,test,testRows),'wrongTarget':metric(wrongs[name],test,testRows),'selection':decision(reg,[features[j] for j in prep[0]],fold['test'],p,test,fold['training'],raw,prep[0])}
   result['noChange']=metric(test['context'],test,testRows);save_file({'observed':test['observed'],'context':test['context']},str(final/'revealed-development-observations.safetensors'));write(dest/'result.json',result);results.append(result);print(fold['id'],rep,{k:v['metrics']['equalGroupRMSE'] for k,v in result['models'].items()},flush=True)
 write(out/'results.json',{'separateTasks':True,'results':results,'biologicalPromotion':False,'nextIndependentEvaluation':'not opened; technical transfer does not establish biological independence'})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('action',choices=['register','prepare','run']);p.add_argument('--inputs',type=Path);p.add_argument('--embeddings',type=Path);p.add_argument('--targets',type=Path);p.add_argument('--binary',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.action=='register':register(a.inputs,a.embeddings,a.targets,a.binary,a.output)
 elif a.action=='prepare':prepare(a.output)
 else:execute(a.output)
