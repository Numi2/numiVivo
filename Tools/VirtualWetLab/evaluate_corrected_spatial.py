#!/usr/bin/env python3
"""Compact exposed-development evaluation of corrected native spatial artifacts."""
import argparse
from pathlib import Path
import numpy as np
import anndata as ad
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,inventory,timestamp
from prepare_spatial_learning import logcounts
from evaluate_spatial_response import metrics

def run(root,kind):
    selection=read(root/'development-selection.json');reg=read(root/'development-registration.json');source=Path(reg['source']);require(not (root/'comparison.json').exists(),'Retain existing evaluation')
    if kind=='receiver':
        old=load_file(str(source/'development-results.safetensors'));pred={v:load_file(str(root/v/'development-prediction/prediction.safetensors'))['mean'] for v in selection['models']};mask=load_file(str(root/'receiver-pretrained/chip2.safetensors'))['mask'][0].astype(bool)&load_file(str(root/'receiver-pretrained/chip3.safetensors'))['mask'][0].astype(bool)
        previous=read(source/'comparison.json');features=read(source/'scored-features.json');write(root/'scored-features.json',features);axis=read(root.parent/'receiver-feature-axis.json');write(root/'features.json',axis);mask=mask&np.isin(axis,features)
        # All these coordinates were already scored in the retained exposed source campaign.
        comparison={'scope':'exposed corrected-optimizer receiving-population development','groups':previous['groups'],'models':{},'biologicalPromotion':False,'spatialInformationBenefit':'UNAVAILABLE: sections unresolved; no verified spatial edges or preservation controls'}
        for name,p in {**pred,'no-change':old['control']}.items():
            err=np.sqrt(np.mean((p[:,mask]-old['observed'][:,mask])**2,axis=1));comparison['models'][name]={'groupRMSE':err.tolist(),'meanGroupRMSE':float(err.mean())}
        save_file({**pred,**{k:old[k] for k in ['control','observed','observedVariance']}},str(root/'development-results.safetensors'))
        registration=read(source/'registration.json');registration.update(createdAt=timestamp(),correctedOptimizer=True,developmentRegistrationSHA256=sha(root/'development-registration.json'),modelVersion='v04-corrected-receiver-adam',binarySHA256=reg['binarySHA256']);write(root/'registration.json',registration)
    else:
        prepared=read(root/'prepared.json');cohort=Path(prepared['cohort']);a=ad.read_h5ad(cohort/'chip1.h5ad');rows=read(root/'chip1-rows.json');query=load_file(str(root/'neighborhood/chip1.safetensors'));mask=a.var.measured_in_source.to_numpy()&load_file(str(root/'neighborhood/chip2.safetensors'))['mask'][0].astype(bool)
        pred={v:load_file(str(root/v/'test-prediction/prediction.safetensors')) for v in selection['models']};baseline=load_file(str(root/'baselines.safetensors'));types=[tuple(x) for x in read(root/'baseline-types.json')]
        pred['no-change']={'mean':query['context']};pred['training-mean']={'mean':query['context']+baseline['prior'][query['target']]}
        pred['cell-type-matched']={'mean':np.asarray([query['context'][i]+(baseline['typePrior'][types.index((r['target'],r['role'],r['cellType']))] if (r['target'],r['role'],r['cellType']) in types else baseline['prior'][r['targetIndex']]) for i,r in enumerate(rows)])}
        comparison={'scope':'exposed corrected-optimizer spatial development','modelSelectionSHA256':sha(root/'development-selection.json'),'groups':[],'aggregate':{},'biologicalPromotion':False,'spatialBenefitEstablished':False,'limits':['One exposed mouse; four isolated controls','Unverified section identities; neighborhood context is a hypothesis','Fixed endpoint; no cell-specific trajectories or propagation']}
        for key in sorted({(r['target'],r['cellType'],r['neighborhood'],r['role']) for r in rows}):
            ix=[i for i,r in enumerate(rows) if (r['target'],r['cellType'],r['neighborhood'],r['role'])==key];ids=sorted({j for i in ix for j in rows[i]['outcomeRows']});obs=logcounts(a,ids);control=query['context'][ix].mean(0);scores={}
            for name,p in pred.items():
                mean=p['mean'][ix].mean(0);variance=np.maximum((p['variance'][ix]+p['mean'][ix]**2).mean(0)-mean**2,1e-4) if 'variance' in p else None
                scores[name]=metrics(mean,variance,obs,control,prepared['metricPanelIndices'],mask)
            comparison['groups'].append({'target':key[0],'cellType':key[1],'neighborhood':key[2],'role':key[3],'cells':len(ids),'metrics':scores})
        for name in pred:
            by={t:float(np.mean([r['metrics'][name]['rmse'] for r in comparison['groups'] if r['target']==t])) for t in sorted({r['target'] for r in rows})};comparison['aggregate'][name]={'equalTargetRMSE':float(np.mean(list(by.values()))),'byTarget':by}
    write(root/'comparison.json',comparison);write(root/'artifact-seal.json',{'createdAt':timestamp(),'files':inventory(root),'binarySHA256':reg['binarySHA256'],'biologicalPromotion':False});print(kind,comparison.get('aggregate',comparison.get('models')),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--kind',choices=['spatial','receiver'],required=True);a=p.parse_args();run(a.inputs,a.kind)
