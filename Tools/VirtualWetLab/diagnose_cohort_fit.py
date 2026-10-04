#!/usr/bin/env python3
"""Registered follow-up: optimization versus descriptor and bag-noise limitations."""
import argparse,time
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,timestamp
from train_cohort_recovery import STEPS,key,measure,predict
from train_intervention_design import native

def run(campaign,out):
    require(not out.exists(),'Retain previous follow-up');out.mkdir();reg=read(campaign/'registration.json');require((campaign/'selection.json').exists(),'Finish original registered budget first');inputs=Path(reg['inputs']);binary=Path(reg['binary'])
    arms=['target-ID','descriptor-high-lr','group-mean-descriptor']
    write(out/'registration.json',{'createdAt':timestamp(),'parentRegistrationSHA256':sha(campaign/'registration.json'),'parentSelectionSHA256':sha(campaign/'selection.json'),'ownerSHA256':sha(__file__),'binarySHA256':sha(binary),'scope':'additional exposed development diagnosis; no new biological evaluation','reason':'Single-study K562 fits well but 64-group Jurkat and full cohort underfit. Distinguish descriptor optimization, target lookup, and bag noise at fixed architecture.','arms':arms,'steps':STEPS,'target-ID':'same network; known mask one, biological target descriptor zero, context metadata retained; not unseen-target qualification','descriptor-high-lr':'one higher constant Adam learning rate .003 versus preregistered .001','group-mean-descriptor':'same .001 model on source/context/target means, masks intersected, observational variance retained; qualification against original bags/group means','selection':'diagnostic only; no replacement of parent validation-selected artifact','stopping':'three arms, 5760 updates each'})
    raw={s:load_file(str(inputs/(s+'.safetensors'))) for s in ('training','validation')};rows={s:read(inputs/(s+'-rows.json')) for s in raw};results={}
    for arm in arms:
        f=out/arm;f.mkdir();plan=read(inputs/'plan.json');plan.update(optimizer='adam',objective='mean',steps=STEPS,trainingBudget=STEPS[-1],sampling='group',diagnostics=False,learningRate=.003 if arm=='descriptor-high-lr' else .001);write(f/'plan.json',plan);data={}
        for split,a in raw.items():
            b={k:v.copy() for k,v in a.items()}
            if arm=='target-ID':b['descriptor'][:,:149]=0;b['known'][:]=1
            save_file(b,str(f/(split+'.safetensors')));data[split]=b
        training=f/'training.safetensors'
        if arm=='group-mean-descriptor':
            a=data['training'];groups=[[i for i,r in enumerate(rows['training']) if key(r)==k] for k in sorted(set(map(key,rows['training'])))];b={}
            for k,v in a.items():
                if k in ('prior','responseScale'):b[k]=v.copy()
                elif k=='stratum':b[k]=np.arange(len(groups),dtype=np.int32)
                elif k in ('target','known','descriptor','contextMask','mask'):b[k]=np.asarray([v[g].min(0) if k in ('mask','contextMask') else v[g[0]] for g in groups],dtype=v.dtype)
                else:b[k]=np.asarray([v[g].mean(0) for g in groups],dtype=v.dtype)
            training=f/'group-means.safetensors';save_file(b,str(training))
        native(binary,'train',f/'plan.json',training,f/'training');checkpoints=[]
        for step in STEPS:
            c={'step':step}
            for split,a in data.items():
                pred,digest=predict(binary,f/'plan.json',f/(split+'.safetensors'),f/'training'/f'weights-{step}.safetensors');c[split]=measure(pred['mean'],a,rows[split]);c[split]['predictionSHA256']=digest
            write(f/f'checkpoint-{step}.json',c);checkpoints.append(c)
        write(f/'result.json',{'checkpoints':checkpoints});results[arm]=checkpoints[-1];print(arm,results[arm]['training']['equalGroupRMSE'],flush=True)
    write(out/'summary.json',results)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.campaign,a.output)
