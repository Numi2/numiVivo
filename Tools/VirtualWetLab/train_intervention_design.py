#!/usr/bin/env python3
"""Execute and seal a fixed multi-study native campaign before outcome reveal."""
import argparse, subprocess
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from wetlab import read,write,sha,require,timestamp,inventory
from design_objective import calibrated_controls

VARIANTS=('target-descriptor','target-ID','shuffled-target-descriptor')

def native(binary,mode,plan,data,out,weights=None):
    command=[str(binary),'spatial-response',mode,str(plan),str(data),str(out)]+([str(weights)] if weights else [])
    with out.with_suffix('.log').open('w') as log:
        code=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT).returncode
    require(code==0,'Native execution failed: '+str(out));return command

def validation_score(mean,observed,rows):
    scores=np.sqrt(((mean-observed)**2).mean(1));by={}
    for value,row in zip(scores,rows):by.setdefault((row['source'],row['target']),[]).append(float(value))
    source={}
    for (s,t),v in by.items():source.setdefault(s,[]).append(float(np.mean(v)))
    return float(np.mean([np.mean(x) for x in source.values()]))

def train(inputs,binary,out):
    require(not out.exists(),'Campaign already started; retain and inspect instead of overwriting');out.mkdir()
    prepared=read(inputs/'prepared.json')
    for f,h in prepared['files'].items():require(sha(inputs/f)==h,'Changed preparation '+f)
    raw={s:load_file(str(inputs/(s+'.safetensors'))) for s in ('training','validation','reserved')}
    meta={s:read(inputs/(s+'-rows.json')) for s in raw};plan=read(inputs/'plan.json');targetmeta=read(inputs/'targets.json');targets=sorted(targetmeta)
    write(out/'registration.json',{'createdAt':timestamp(),'inputs':str(inputs.resolve()),'preparedSHA256':sha(inputs/'prepared.json'),
      'binary':str(binary.resolve()),'binarySHA256':sha(binary),'ownerSHA256':sha(__file__),'planSHA256':sha(inputs/'protocol.json'),
      'variants':list(VARIANTS),'heldOutResponseRead':False,'pretrainedExpressionWeights':False})
    rng=np.random.default_rng(271828);permutation=rng.permutation(len(targets));descs=np.asarray([targetmeta[t]['descriptor'] for t in targets],np.float32)
    known=set(raw['training']['target'].tolist());delta=raw['training']['observed']-raw['training']['context']
    baseline=np.zeros((len(targets),plan['featureCount']),np.float32)
    for i in known:baseline[i]=delta[raw['training']['target']==i].mean(0)
    # Target-agnostic training mean is a fixed molecular ranking baseline. ID
    # lookup receives the identical native network and no empirical prior.
    baseline_query=np.asarray([delta[np.asarray([r['modality']==q['modality'] for r in meta['training']])].mean(0) for q in meta['reserved']],np.float32)
    predictions={'no-change':raw['reserved']['context'].copy(),'training-mean':raw['reserved']['context']+baseline_query}
    if read(inputs/'protocol.json').get('targetMatchedBaseline'):
        matched=[]
        for q,fallback in zip(meta['reserved'],baseline_query):
            selected=np.asarray([r['target']==q['target'] and r['modality']==q['modality'] and r['context']==q['context'] for r in meta['training']])
            matched.append(delta[selected].mean(0) if selected.any() else fallback)
        predictions['training-mean']=raw['reserved']['context']+np.asarray(matched,np.float32)
    selections={}
    for variant in VARIANTS:
        folder=out/variant;folder.mkdir();write(folder/'plan.json',plan)
        for split,a in raw.items():
            data={k:v.copy() for k,v in a.items()}
            if variant=='target-ID':
                data['descriptor'][:,:149]=0;data['known'][:,0]=np.isin(data['target'],list(known)).astype(np.float32)
            elif variant=='shuffled-target-descriptor':data['descriptor'][:,:149]=descs[permutation[data['target']]]
            save_file(data,str(folder/(split+'.safetensors')))
        native(binary,'train',folder/'plan.json',folder/'training.safetensors',folder/'training')
        scores=[]
        for step in plan['steps']:
            dest=folder/f'validation-{step}';native(binary,'predict',folder/'plan.json',folder/'validation.safetensors',dest,folder/'training'/f'weights-{step}.safetensors')
            p=load_file(str(dest/'prediction.safetensors'))
            scores.append({'step':step,'score':validation_score(p['mean'],raw['validation']['observed'],meta['validation'])})
        selected=min(scores,key=lambda x:(x['score'],x['step']));selections[variant]={'scores':scores,'selected':selected}
        write(folder/'selection.json',selections[variant]);weights=folder/'training'/f"weights-{selected['step']}.safetensors"
        native(binary,'predict',folder/'plan.json',folder/'reserved.safetensors',folder/'reserved-prediction',weights)
        predictions[variant]=load_file(str(folder/'reserved-prediction/prediction.safetensors'))['mean']
        print(variant,selected,flush=True)
    # Same biological descriptors, simple regularized response comparator.
    pca=PCA(n_components=min(16,len(raw['training']['context'])-1),svd_solver='full').fit(raw['training']['context'])
    tx=np.c_[raw['training']['descriptor'],pca.transform(raw['training']['context'])];qx=np.c_[raw['reserved']['descriptor'],pca.transform(raw['reserved']['context'])]
    center=tx.mean(0);scale=np.maximum(tx.std(0),.1);model=Ridge(alpha=10).fit((tx-center)/scale,delta)
    predictions['regularized-descriptor-ridge']=raw['reserved']['context']+model.predict((qx-center)/scale)
    save_file({'coefficient':model.coef_.astype(np.float32),'intercept':model.intercept_.astype(np.float32),'center':center.astype(np.float32),'scale':scale.astype(np.float32),'pcaMean':pca.mean_.astype(np.float32),'pcaComponents':pca.components_.astype(np.float32)},str(out/'ridge-weights.safetensors'))
    # Transparent GO-neighbor response transfer: reproduction of the relationship
    # principle, NOT the trained GEARS neural network or a published GEARS score.
    transfer=[];neighbors=[]
    for q in meta['reserved']:
        possible=sorted({r['target'] for r in meta['training'] if r['modality']==q['modality']});go=set(targetmeta[q['target']]['go']);rank=[]
        for target in possible:
            other=set(targetmeta[target]['go']);rank.append((len(go&other)/max(1,len(go|other)),target))
        rank=sorted(rank,key=lambda x:(-x[0],x[1]))[:3];weights=np.asarray([x[0] for x in rank]);weights=weights/weights.sum() if weights.sum()>0 else np.ones(len(rank))/len(rank)
        values=[delta[[r['target']==t and r['modality']==q['modality'] for r in meta['training']]].mean(0) for _,t in rank]
        transfer.append(np.average(values,axis=0,weights=weights));neighbors.append({'target':q['target'],'source':q['source'],'neighbors':rank})
    predictions['GO-neighbor-transfer']=raw['reserved']['context']+np.asarray(transfer,np.float32)
    write(out/'GO-transfer-neighbors.json',neighbors);save_file({k:v.astype(np.float32) for k,v in predictions.items()},str(out/'all-predictions.safetensors'))
    # Calibration data are training-only, reduced to one mean per source-target.
    means=[]
    for key in sorted({(r['source'],r['target']) for r in meta['training']}):means.append(delta[[ (r['source'],r['target'])==key for r in meta['training']]].mean(0))
    write(out/'metric-calibration.json',calibrated_controls(means))
    write(out/'prediction-seal.json',{'createdAt':timestamp(),'allCandidatesSealedTogether':True,'files':inventory(out),
      'queriesSHA256':sha(inputs/'reserved.safetensors'),'metadataSHA256':sha(inputs/'reserved-rows.json'),
      'protocolSHA256':sha(inputs/'protocol.json'),'unseenTargets':sorted(set(r['target'] for r in meta['reserved'])-set(r['target'] for r in meta['training']))})
    print('SEALED; reserved responses not opened',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();train(a.inputs,a.binary,a.output)
