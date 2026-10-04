#!/usr/bin/env python3
"""Bounded fixed-architecture DEVELOPMENT campaign. Never opens reserved outcomes.

All checkpoint weights, row/group metrics and exact inputs are retained. Prediction
scratch tensors are reduced to diagnostics and hashed, then removed; replay uses
retained weights and inputs. Validation is a technical capture from ONE study.
"""
import argparse, json, tempfile, time
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file, save_file
from wetlab import read, write, sha, require, timestamp
from train_intervention_design import native

STEPS=[240,720,1440,2880,5760]

def key(r): return (r['source'],r['context'],r['target'])
def measure(mean,a,rows):
    records=[];pred=[];truth=[]
    for k in sorted(set(map(key,rows))):
        ix=[i for i,r in enumerate(rows) if key(r)==k]; mask=a['mask'][ix].min(0)>0
        y=(a['observed'][ix]-a['context'][ix]).mean(0);p=(mean[ix]-a['context'][ix]).mean(0)
        supported=mask & (np.abs(y)>=.1)
        records.append({'source':k[0],'condition':k[1],'target':k[2],'bags':len(ix),'features':int(mask.sum()),
          'responseRMSE':float(np.sqrt(np.mean((p[mask]-y[mask])**2))),
          'noChangeRMSE':float(np.sqrt(np.mean(y[mask]**2))),
          'directionAccuracy':float(np.mean(np.sign(p[supported])==np.sign(y[supported]))) if supported.any() else None,
          'directionFeatures':int(supported.sum()),'predictedResponseRMS':float(np.sqrt(np.mean(p[mask]**2)))})
        pred.append(p);truth.append(y)
    # Identity discrimination is within source and condition, never across assay scales.
    for i,r in enumerate(records):
        candidates=[j for j,s in enumerate(records) if (s['source'],s['condition'])==(r['source'],r['condition'])]
        mask=a['mask'][[j for j,x in enumerate(rows) if key(x)==(r['source'],r['condition'],r['target'])][0]]>0
        best=min(candidates,key=lambda j:float(np.mean((pred[i][mask]-truth[j][mask])**2)))
        r['nearestMeasuredTarget']=records[best]['target'];r['identityCorrect']=best==i
    bySource={s:float(np.mean([r['responseRMSE'] for r in records if r['source']==s])) for s in sorted(set(r['source'] for r in records))}
    return {'groups':records,'bySource':bySource,'equalSourceRMSE':float(np.mean(list(bySource.values()))),
      'equalGroupRMSE':float(np.mean([r['responseRMSE'] for r in records])),
      'identityAccuracy':float(np.mean([r['identityCorrect'] for r in records])),
      'directionAccuracy':float(np.mean([r['directionAccuracy'] for r in records if r['directionAccuracy'] is not None])),
      'predictedBetweenGroupSD':float(np.std(pred,axis=0).mean()),'observedBetweenGroupSD':float(np.std(truth,axis=0).mean())}

def predict(binary,plan,data,weights):
    with tempfile.TemporaryDirectory(prefix='cohort-predict-') as tmp:
        out=Path(tmp)/'prediction';native(binary,'predict',plan,data,out,weights)
        return load_file(str(out/'prediction.safetensors')),sha(out/'prediction.safetensors')

def register(inputs,binary,out):
    require(not out.exists(),'Registration exists');out.mkdir()
    rows=read(inputs/'training-rows.json');sources=sorted(set(r['source'] for r in rows));arms=[]
    for source in sources:
        for objective in ('mean','distribution'):
            arms.append(dict(id=source+'-'+objective,source=source,objective=objective,sampling='group',learningRate=.001,condition='metadata'))
    for objective in ('mean','distribution'):
        for sampling in ('group','source'):
            arms.append(dict(id='all-'+objective+'-'+sampling,source=None,objective=objective,sampling=sampling,learningRate=.001,condition='metadata'))
    arms += [dict(id='all-mean-source-low-lr',source=None,objective='mean',sampling='source',learningRate=.0003,condition='metadata')]
    for condition in ('RNA-only','shuffled-metadata'):
        arms.append(dict(id='all-mean-source-'+condition,source=None,objective='mean',sampling='source',learningRate=.001,condition=condition))
    write(out/'registration.json',{'createdAt':timestamp(),'scope':'exposed development only; no independent validation or new reservation',
      'inputs':str(inputs.resolve()),'inputsSHA256':{f:sha(inputs/f) for f in ['plan.json','training.safetensors','training-rows.json','validation.safetensors','validation-rows.json']},
      'binary':str(binary.resolve()),'binarySHA256':sha(binary),'ownerSHA256':sha(__file__),'architecture':'unchanged width 64 cell/descriptor/target fusion',
      'priorQualification':'v04-recovery diagnosis-v2: qualified three-target Adam training, not biological validation',
      'steps':STEPS,'trainingBudget':STEPS[-1],'arms':arms,'groups':len(set(map(key,rows))),
      'selection':'all-source metadata arms only: minimum technical-validation equal-source RMSE, ties earlier checkpoint then arm ID; report training fit at selected and final checkpoints',
      'fitCriteria':'Lower group response RMSE than no change; target identity and direction reported jointly. Response variation is diagnostic only, never a selection criterion.',
      'directionThreshold':.1,'validationLimit':'Only GSE92872 stimulated technical capture6; no validation available for GSE90546 or GSE90063, not independent animals',
      'conditionProbes':'Development GSE92872 stimulated/unstimulated paired targets. Swap control RNA only, metadata only, and both; correct paired condition outcomes are already exposed development. No GSM2406677 use.',
      'uncertainty':'Descriptive group summaries only; bags/genes are not biological replicates','stopping':'Exactly 13 arms, 5760 updates each, no architecture or feature changes; retain failed arms'})

def execute(out):
    reg=read(out/'registration.json');inputs=Path(reg['inputs']);binary=Path(reg['binary'])
    require(sha(binary)==reg['binarySHA256'] and sha(__file__)==reg['ownerSHA256'],'Changed registered executable/owner')
    for f,h in reg['inputsSHA256'].items():require(sha(inputs/f)==h,'Changed input '+f)
    raw={s:load_file(str(inputs/(s+'.safetensors'))) for s in ('training','validation')};metadata={s:read(inputs/(s+'-rows.json')) for s in raw};results={}
    for arm in reg['arms']:
        folder=out/arm['id']
        if (folder/'result.json').exists():results[arm['id']]=read(folder/'result.json');continue
        require(not folder.exists(),'Interrupted arm must be retained; no silent restart: '+str(folder));folder.mkdir()
        plan=read(inputs/'plan.json');plan.update({k:arm[k] for k in ('objective','sampling','learningRate')});plan.update(steps=STEPS,trainingBudget=STEPS[-1],diagnostics=False)
        write(folder/'plan.json',plan);data={};rows={}
        for split,a in raw.items():
            ix=np.asarray([i for i,r in enumerate(metadata[split]) if not arm['source'] or r['source'].split('-capture')[0]==arm['source']],int)
            if not len(ix):continue
            b={k:(v.copy() if k in ('prior','responseScale') else v[ix].copy()) for k,v in a.items()};rr=[metadata[split][i] for i in ix]
            sourceNames=sorted(set(r['source'] for r in rr));b['source']=np.asarray([sourceNames.index(r['source']) for r in rr],np.int32)
            if arm['condition']=='RNA-only':b['descriptor'][:,-2:]=0
            elif arm['condition']=='shuffled-metadata':
                # Only shuffle supported condition labels WITHIN source, preserving modality/species.
                rng=np.random.default_rng(271828)
                for s in sourceNames:
                    ids=np.flatnonzero(b['source']==sourceNames.index(s));b['descriptor'][ids,-2:]=b['descriptor'][rng.permutation(ids),-2:]
            save_file(b,str(folder/(split+'.safetensors')));write(folder/(split+'-rows.json'),[{k:r[k] for k in ('source','context','target','unit')} for r in rr]);data[split]=b;rows[split]=rr
        start=time.monotonic();native(binary,'train',folder/'plan.json',folder/'training.safetensors',folder/'training');checkpoints=[]
        for step in STEPS:
            result={'step':step,'weightsSHA256':sha(folder/'training'/f'weights-{step}.safetensors')}
            for split,a in data.items():
                pred,digest=predict(binary,folder/'plan.json',folder/(split+'.safetensors'),folder/'training'/f'weights-{step}.safetensors')
                result[split]=measure(pred['mean'],a,rows[split]);result[split]['predictionSHA256']=digest
            result.setdefault('validation',{'status':'UNAVAILABLE: no withheld development capture for this source'})
            checkpoints.append(result);write(folder/f'checkpoint-{step}.json',result)
        result={'arm':arm,'seconds':time.monotonic()-start,'checkpoints':checkpoints};write(folder/'result.json',result);results[arm['id']]=result
        print(arm['id'],checkpoints[-1]['training']['equalGroupRMSE'],flush=True)
    candidates=[(c['validation']['equalSourceRMSE'],c['step'],name) for name,r in results.items() if r['arm']['source'] is None and r['arm']['condition']=='metadata' for c in r['checkpoints']]
    score,step,name=min(candidates);write(out/'selection.json',{'arm':name,'step':step,'validationRMSE':score,'biologicalPromotion':False,'modelSelectionFrozenAt':timestamp(),'validationScope':reg['validationLimit']})
    write(out/'summary.json',{k:{'final':r['checkpoints'][-1],'seconds':r['seconds']} for k,r in results.items()})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['register','run']);p.add_argument('--inputs',type=Path);p.add_argument('--binary',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.action=='register':register(a.inputs,a.binary,a.output)
    else:execute(a.output)
