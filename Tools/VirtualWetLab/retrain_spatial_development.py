#!/usr/bin/env python3
"""Corrected optimizer on retained spatial/receiver inputs; development only."""
import argparse, subprocess
from pathlib import Path
from safetensors.numpy import load_file
from wetlab import read,write,sha,require,timestamp
from train_intervention_design import native
from train_cohort_recovery import measure
from bounded_spatial_prediction import predict

def clone(a,b):subprocess.run(['cp','-c',str(a),str(b)],check=True)
def run(source,binary,out,kind,resume=False):
    if resume:
        require(out.exists(),'No retained attempt');registered=read(out/'development-registration.json');require(registered['binarySHA256']==sha(binary),'Runtime changed')
        for f,h in registered['inputsSHA256'].items():require(sha(source/f)==h,'Input changed')
        write(out/('scoring-recovery-'+str(len(list(out.glob('scoring-recovery*.json')))+1)+'.json'),{'createdAt':timestamp(),'previousOwnerSHA256':registered['ownerSHA256'],'ownerSHA256':sha(__file__),'scoringOwnerSHA256':sha(Path(__file__).with_name('bounded_spatial_prediction.py')),'reason':'Resume stopped scoring/input interruption; retain completed weights and selections, use batch-aligned scoring. See source-restoration and original failure records.','retainedWeights':{str(p.relative_to(out)):sha(p) for p in out.glob('*/training/weights-*.safetensors')},'trainingBudgetUnchanged':True,'predictionHashFormat':'sha256 sorted tensor names then concatenated float32 bytes'})
    else:
        require(not out.exists(),'Retain previous attempt');out.mkdir()
    variants=['neighborhood','no-neighborhood','shuffled-neighborhood'] if kind=='spatial' else ['receiver-pretrained','anchor-reference-pretrained','receiver-spatial-only']
    steps=[240,1440,4320]
    if not resume:write(out/'development-registration.json',{'createdAt':timestamp(),'source':str(source),'sourceRegistrationSHA256':sha(source/('preregistration.json' if kind=='spatial' else 'registration.json')),'binarySHA256':sha(binary),'ownerSHA256':sha(__file__),'inputsSHA256':{str(Path(v)/f):sha(source/v/f) for v in variants for f in ['plan.json','chip1.safetensors','chip2.safetensors','chip3.safetensors']},'scope':'permanently exposed spatial development, not new validation','variants':variants,'steps':steps,'optimizer':'adam','objective':'distribution','architecture':'unchanged','selection':'minimum chip3 equal-target RMSE; earlier step for ties','edges':'unverified sections; spatial descriptors remain hypotheses, preservation blocked','stopping':'three variants, 4320 updates each'})
    if kind=='spatial':
        for name in ['prepared.json','preregistration.json','features.json','chip1-rows.json','chip2-rows.json','chip3-rows.json','baselines.safetensors','baseline-types.json']:
            if not (out/name).exists():clone(source/name,out/name)
    else:
        for name in ['chip1-rows.json','chip2-rows.json','chip3-rows.json','warm-start.safetensors']:
            if not (out/name).exists():clone(source/name,out/name)
    results={}
    for variant in variants:
        completed=out/variant
        if resume and (completed/'selection.json').exists() and (completed/('test-prediction' if kind=='spatial' else 'development-prediction')/'prediction.safetensors').exists():
            results[variant]=read(completed/'selection.json')['selected'];continue
        folder=out/variant;folder.mkdir(exist_ok=resume);old=source/variant
        plan=read(old/'plan.json');plan.update(steps=steps,trainingBudget=steps[-1],optimizer='adam',objective='distribution',diagnostics=False)
        if not (folder/'plan.json').exists():write(folder/'plan.json',plan)
        else:require(read(folder/'plan.json')==plan,'Plan changed during recovery')
        for chip in ('chip1','chip2','chip3'):
            if not (folder/(chip+'.safetensors')).exists():clone(old/(chip+'.safetensors'),folder/(chip+'.safetensors'))
        warm=out/'warm-start.safetensors' if kind=='receiver' and variant!='receiver-spatial-only' else None
        if not (folder/'training'/f'weights-{steps[-1]}.safetensors').exists():native(binary,'train',folder/'plan.json',folder/'chip2.safetensors',folder/'training',warm)
        checkpoints=[]
        for step in steps:
            result={'step':step,'weightsSHA256':sha(folder/'training'/f'weights-{step}.safetensors')}
            for split,chip in [('training','chip2'),('validation','chip3')]:
                a=load_file(str(folder/(chip+'.safetensors')));rows=read(out/(chip+'-rows.json'))
                rr=[{'source':chip,'context':r.get('cellType',r.get('receivingType','unknown'))+'|'+r['role'],'target':r['target']} for r in rows]
                prediction,digest=predict(binary,folder/'plan.json',folder/(chip+'.safetensors'),folder/'training'/f'weights-{step}.safetensors')
                result[split]=measure(prediction['mean'],a,rr);result[split]['predictionSHA256']=digest
                # Match the original equal-target selection, without substituting a new metric.
                import numpy as np
                by=[]
                for target in sorted(set(r['target'] for r in rows)):
                    ix=[i for i,r in enumerate(rows) if r['target']==target];mask=a['mask'][ix];by.append(float(np.sqrt((((prediction['mean'][ix]-a['observed'][ix])**2)*mask).sum()/mask.sum())))
                result[split]['equalTargetRMSE']=float(np.mean(by))
            checkpoints.append(result);write(folder/f'checkpoint-{step}.json',result)
        selected=min(checkpoints,key=lambda x:(x['validation']['equalTargetRMSE'],x['step']))
        choice={'step':selected['step'],'equalTargetRMSE':selected['validation']['equalTargetRMSE']}
        write(folder/'selection.json',{'selected':choice,'candidates':checkpoints,'biologicalPromotion':False,'scope':'exposed development'})
        dest=folder/('test-prediction' if kind=='spatial' else 'development-prediction');native(binary,'predict',folder/'plan.json',folder/'chip1.safetensors',dest,folder/'training'/f"weights-{choice['step']}.safetensors")
        results[variant]=choice;print(kind,variant,choice,flush=True)
    write(out/'development-selection.json',{'models':results,'modelSelectionFrozenAt':timestamp(),'biologicalPromotion':False})
    if kind=='spatial':
        from learned_spatial import build
        config=build(out,binary,out/'assay');a=read(config);a.update(id='spatial-response-v04-corrected',title='v0.4 corrected Adam · exposed spatial development',modelVersion='v04-cohort-adam',conditions=[{'id':'measured-endpoint','title':'Source study measured endpoint; dose/time manipulation unavailable'}]);a['limits'].append('Corrected optimizer development model; pooled neighbor target is not receiver-specific biology. No validated spatial benefit.')
        # New artifact identity; never modify an existing published assay.
        config.write_text(__import__('json').dumps(a,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--kind',choices=['spatial','receiver'],required=True);p.add_argument('--resume',action='store_true');a=p.parse_args();run(a.source,a.binary,a.output,a.kind,a.resume)
