#!/usr/bin/env python3
"""Bounded development-only target-learning diagnosis on the existing native model."""
import argparse
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from sklearn.linear_model import Ridge
from wetlab import read,write,sha,require,inventory,timestamp
from train_intervention_design import native

ARMS=[('legacy','distribution','sgd',.0001,False,False),('mean-sgd-regularized','mean','sgd',.0001,False,False),('mean-sgd','mean','sgd',0.,False,False),('mean-adam','mean','adam',0.,False,False),('mean-adam-scaled','mean','adam',0.,True,False),('mean-adam-ID','mean','adam',0.,False,True),('mean-adam-regularized','mean','adam',.0001,False,False),('distribution-adam','distribution','adam',.0001,False,False)]

def prepare(inputs,out):
    require(not out.exists(),'Retain prior diagnosis');out.mkdir();raw=load_file(str(inputs/'training.safetensors'));rows=read(inputs/'training-rows.json')
    keys=sorted({r['target'] for r in rows if r['source']=='GSE92872-training' and r['context']=='stimulated'})
    if not keys:keys=sorted({r['target'] for r in rows if r['source'].startswith('GSE92872') and r['context']=='stimulated'})
    groups=[]
    for t in keys:
        ix=[i for i,r in enumerate(rows) if r['target']==t and r['source'].startswith('GSE92872') and r['context']=='stimulated'];r=rows[ix[0]]
        if r['outcomeCellCount']>=100 and r['controlCellCount']>=100:groups.append((r['outcomeCellCount'],t,ix))
    chosen=sorted(groups,key=lambda x:(-x[0],x[1]))[:3];require(len(chosen)==3,'Need three supported development targets');a={k:[] for k in ['context','descriptor','target','known','observed','observedVariance','mask','stratum']}
    for j,(_,t,ix) in enumerate(chosen):
        for k in a:a[k].append(j if k=='stratum' else raw[k][ix[0]] if k in ('context','descriptor','target','known','mask') else raw[k][ix].mean(0))
    a={k:np.asarray(v,np.int32 if k in ('target','stratum') else np.float32) for k,v in a.items()};a['prior']=raw['prior'];save_file(a,str(out/'training.safetensors'));write(out/'rows.json',[rows[x[2][0]] for x in chosen]);write(out/'features.json',read(inputs/'features.json'))
    validation=load_file(str(inputs/'validation.safetensors'));vr=read(inputs/'validation-rows.json');ix=[i for i,r in enumerate(vr) if r['context']=='stimulated' and r['target'] in {t for _,t,_ in chosen}];save_file({k:v[ix] for k,v in validation.items()},str(out/'validation.safetensors'));write(out/'validation-rows.json',[vr[i] for i in ix]);write(out/'base-plan.json',read(inputs/'plan.json'))
    write(out/'preregistration.json',{'createdAt':timestamp(),'scope':'development-only engineering qualification; all observations previously exposed','sourceInputs':str(inputs),'preparedSHA256':sha(inputs/'prepared.json'),'groups':'three largest source-target groups with >=100 assigned cells and >=100 explicit controls in GSE92872 stimulated training captures; targets selected by counts, not effects','architecture':'unchanged width64 VivoCellResponseMLXModel; no new layers','targets':[t for _,t,_ in chosen],'arms':ARMS,'budget':'one seed271828; eight arms; checkpoints240,720,1440; no retries or size search','fitGate':'training RMSE <=20% of no-change AND correct nearest response identity for all three groups AND swapped target inputs reproduce corresponding predictions to maxAbs<=1e-5','scaling':'per-feature training response SD, floor0.1; never validation statistics','gradientQualification':'finite differences on maximal decoder bias coordinate, masked observation corruption invariant, per-layer gradient/update norms at steps1,2,10,1440','regressionComparator':'descriptor-conditioned ridge alpha1e-6 on identical three groups; engineering comparator, not native biological promotion','generalization':'only after fit: exposed capture6 reports target-mean and no-change comparison; no biological validation claim','reservation':'no new evaluation cohort reserved or opened until diagnosis resolved','oldFailuresPreserved':True})

def fit(out,binary):
    require(not (out/'results.json').exists(),'Diagnosis already evaluated');raw=load_file(str(out/'training.safetensors'));val=load_file(str(out/'validation.safetensors'));rows=read(out/'rows.json');vr=read(out/'validation-rows.json');base=read(out/'base-plan.json');delta=raw['observed']-raw['context'];noise=float(np.sqrt(np.mean(delta**2)));results={}
    x=raw['descriptor'][:,:149];ridge=Ridge(alpha=1e-6).fit(x,delta);rp=raw['context']+ridge.predict(x)
    results['descriptor-regression']={'trainingRMSE':float(np.sqrt(np.mean((rp-raw['observed'])**2))),'scope':'same-input regularized linear engineering comparator'}
    for name,objective,optimizer,decay,scaled,ident in ARMS:
        f=out/name;f.mkdir();plan={**base,'objective':objective,'optimizer':optimizer,'weightDecay':decay,'diagnostics':True};write(f/'plan.json',plan);a={k:v.copy() for k,v in raw.items()};b={k:v.copy() for k,v in val.items()}
        for q in (a,b):
            if ident:q['descriptor'][:,:149]=0;q['known'][:]=1
            if scaled:q['responseScale']=np.maximum(delta.std(0),.1).astype(np.float32)
        save_file(a,str(f/'training.safetensors'));save_file(b,str(f/'validation.safetensors'));native(binary,'train',f/'plan.json',f/'training.safetensors',f/'training')
        checks=[]
        for step in plan['steps']:
            native(binary,'predict',f/'plan.json',f/'training.safetensors',f/f'fit-{step}',f/'training'/f'weights-{step}.safetensors');p=load_file(str(f/f'fit-{step}/prediction.safetensors'))['mean'];rmse=float(np.sqrt(np.mean((p-raw['observed'])**2)));identity=float(np.mean(np.argmin(((p[:,None]-raw['observed'][None,:])**2).mean(2),axis=1)==np.arange(len(p))));checks.append({'step':step,'trainingRMSE':rmse,'identityAccuracy':identity,'responseSD':float((p-raw['context']).std(0).mean())})
        step=plan['steps'][-1];query={k:v.copy() for k,v in a.items() if k in ('context','descriptor','target','known')};query['descriptor'][:,:149]=np.roll(query['descriptor'][:,:149],1,axis=0);query['target']=np.roll(query['target'],1);save_file(query,str(f/'swapped.safetensors'));native(binary,'predict',f/'plan.json',f/'swapped.safetensors',f/'swapped',f/'training'/f'weights-{step}.safetensors');sw=load_file(str(f/'swapped/prediction.safetensors'))['mean'];permutation_error=float(np.max(np.abs(sw-np.roll(p,1,axis=0))))
        native(binary,'predict',f/'plan.json',f/'validation.safetensors',f/'validation',f/'training'/f'weights-{step}.safetensors');vp=load_file(str(f/'validation/prediction.safetensors'))['mean'];vd=b['observed']-b['context'];matched=np.asarray([delta[next(i for i,r in enumerate(rows) if r['target']==v['target'])] for v in vr]);di=read(f/'training/training-diagnostics.json');require(all(d['maskedObservationInvariant'] for d in di),'Masked observations affected loss');results[name]={'checkpoints':checks,'targetSwapMaxAbs':permutation_error,'noChangeTrainingRMSE':noise,'observedResponseSD':float(delta.std(0).mean()),'fitGate':checks[-1]['trainingRMSE']<=.2*noise and checks[-1]['identityAccuracy']==1 and permutation_error<=1e-5,'developmentGeneralization':{'RMSE':float(np.sqrt(np.mean((vp-b['observed'])**2))),'noChange':float(np.sqrt(np.mean(vd**2))),'matchedTargetMean':float(np.sqrt(np.mean((vd-matched)**2))),'scope':'exposed technical capture6, not biological validation'}};print(name,results[name],flush=True)
    write(out/'results.json',{'results':results,'biologicalPromotion':False,'trainingFitOnly':True,'binarySHA256':sha(binary),'ownerSHA256':sha(__file__)});write(out/'artifact-seal.json',{'files':inventory(out)})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','fit']);p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--binary',type=Path);a=p.parse_args();prepare(a.inputs,a.output) if a.mode=='prepare' else fit(a.output,a.binary)
