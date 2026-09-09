#!/usr/bin/env python3
"""Train-only donor-response baselines from exact real pseudobulk axes."""
import argparse, hashlib, importlib.metadata, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.kernel_ridge import KernelRidge


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def write(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def prepare_training(raw_control,raw_treated,ids):
    # Only training donor matrices can enter this function.
    control=np.log1p(raw_control/raw_control.sum(axis=1,keepdims=True)*1e6)
    treated=np.log1p(raw_treated/raw_treated.sum(axis=1,keepdims=True)*1e6)
    delta=treated-control
    selected=np.flatnonzero(((raw_control+raw_treated).sum(axis=0)>=10)&(((raw_control+raw_treated)>0).sum(axis=0)>=2))
    assert len(selected)>0
    center=control[:,selected].mean(axis=0);scale=control[:,selected].std(axis=0);scale[scale==0]=1
    context=(control[:,selected]-center)/scale/np.sqrt(len(selected))
    mean=delta.mean(axis=0);median=np.median(delta,axis=0)
    reference=KernelRidge(alpha=1,kernel='linear').fit(context,delta-mean)
    dual=np.linalg.solve(context@context.T+np.eye(len(context)),delta-mean)
    np.testing.assert_allclose(dual,reference.dual_coef_,rtol=1e-10,atol=1e-10)
    top=np.array(sorted(range(len(ids)),key=lambda i:(-abs(mean[i]),ids[i]))[:min(200,len(ids))],dtype=np.int64)
    return dict(selected=selected,center=center,scale=scale,context=context,mean=mean,median=median,dual=dual,top=top),reference


def predict(model,raw_control):
    # No treated outcomes, donor labels or global normalization factors are accepted.
    control=np.log1p(raw_control/raw_control.sum()*1e6)
    context=(control[model['selected']]-model['center'])/model['scale']/np.sqrt(len(model['selected']))
    raw={'noChange':np.zeros_like(control),'meanResponse':model['mean'],'medianResponse':model['median'],
         'contextRidge':model['mean']+(context@model['context'].T)@model['dual']}
    return control,context,raw,{k:np.maximum(control+v,0) for k,v in raw.items()}


def metrics(pred,truth):
    residual=pred-truth
    value={'features':len(truth),'rmse':float(np.sqrt(np.mean(residual**2))),'mae':float(np.mean(abs(residual))),
           'pearson':None,'nonzeroTruthFeatures':int(np.count_nonzero(truth)),'signAgreement':None,'explainedSumSquaresVsZero':None}
    if np.std(pred)>0 and np.std(truth)>0:value['pearson']=float(np.corrcoef(pred,truth)[0,1])
    active=truth!=0
    if active.any():value['signAgreement']=float(np.mean(np.sign(pred[active])==np.sign(truth[active])))
    denominator=np.sum(truth**2)
    if denominator>0:value['explainedSumSquaresVsZero']=float(1-np.sum(residual**2)/denominator)
    return value


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--control',required=True);p.add_argument('--treated',required=True);p.add_argument('--study',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    source=pd.read_csv(a.input/'counts.tsv',sep='\t',index_col=0)
    samples=pd.read_csv(a.input/'samples.tsv',sep='\t').set_index('sampleID')
    assert source.index.is_unique and samples.index.is_unique and list(source.columns)==list(samples.index)
    counts=source.to_numpy(dtype=np.float64).T
    assert np.isfinite(counts).all() and (counts>=0).all() and np.array_equal(counts,np.floor(counts))
    np.testing.assert_array_equal(counts.sum(axis=1),samples.libraryCounts.to_numpy())
    assert (counts.sum(axis=1)>0).all() and set(samples.condition)=={a.control,a.treated}
    donors=sorted(set(samples.donor));assert len(donors)>=3
    control_rows=[];treated_rows=[]
    for donor in donors:
        for condition,rows in [(a.control,control_rows),(a.treated,treated_rows)]:
            match=np.flatnonzero((samples.donor.to_numpy()==donor)&(samples.condition.to_numpy()==condition));assert len(match)==1;rows.append(match[0])
    control=counts[control_rows];treated=counts[treated_rows];ids=np.asarray(source.index,dtype=str)
    prior=json.loads((a.input/'input.json').read_text());markers=[x for x in prior['expectedGenes'] if x in source.index]
    write(a.out/'inputs.json',{'study':a.study,'inputHashes':{f:sha(a.input/f) for f in ['counts.tsv','samples.tsv','input.json']},'source':prior,'donors':donors,'features':len(ids),'markers':markers,'packages':{n:importlib.metadata.version(n) for n in ['numpy','pandas','scikit-learn']},'normalization':'per-library natural-log(1+CPM), complete source gene denominator; original DE size factors ignored'})
    results=[]
    for i,donor in enumerate(donors):
        root=a.out/donor;root.mkdir();keep=np.arange(len(donors))!=i
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always');model,reference=prepare_training(control[keep],treated[keep],ids)
        baseline,query_context,raw,predicted=predict(model,control[i])
        np.testing.assert_allclose(raw['contextRidge'],model['mean']+reference.predict(query_context[None,:])[0],rtol=1e-10,atol=1e-10)
        # Change the sealed treated outcome; training and predictions must stay identical.
        changed=treated.copy();changed[i]=changed[i][::-1]+17
        alternate,_=prepare_training(control[keep],changed[keep],ids)
        for key in model:np.testing.assert_array_equal(model[key],alternate[key])
        _,_,_,second=predict(alternate,control[i])
        for key in predicted:np.testing.assert_array_equal(predicted[key],second[key])
        np.savez_compressed(root/'model.npz',**model,featureIDs=ids,trainingDonors=np.asarray(donors,dtype=str)[keep])
        with np.load(root/'model.npz',allow_pickle=False) as frozen:
            _,_,_,reloaded=predict(frozen,control[i])
            for key in predicted:np.testing.assert_array_equal(predicted[key],reloaded[key])
        # The held-out treated library becomes visible only for evaluation below.
        observed=np.log1p(treated[i]/treated[i].sum()*1e6);truth=observed-baseline
        panels={'allGenes':np.arange(len(ids)),'trainingExpressed':model['selected'],'trainingTop200Response':model['top'],'priorMarkers':np.array([source.index.get_loc(m) for m in markers],dtype=np.int64)}
        fold={'heldOutDonor':donor,'trainingDonors':np.asarray(donors)[keep].tolist(),'warnings':[str(w.message) for w in caught],'models':{},'compositionDiagnostics':{}}
        for name,estimate in predicted.items():
            impliedCPM=float(np.expm1(estimate).sum())
            fold['compositionDiagnostics'][name]={'impliedCPMSum':impliedCPM,'relativeDeviationFromMillion':impliedCPM/1e6-1,'renormalized':False}
            fold['models'][name]={panel:{'response':metrics((estimate-baseline)[ix],truth[ix]),'expression':metrics(estimate[ix],observed[ix])} for panel,ix in panels.items() if len(ix)}
            np.savez_compressed(root/(name+'.npz'),featureIDs=ids,control=baseline,observedTreated=observed,unclippedResponse=raw[name],predictedTreated=estimate,predictedResponse=estimate-baseline)
        write(root/'metrics.json',fold);results.append(fold)
        print(json.dumps({'donor':donor,'trainingTop200ResponseRMSE':{k:v['trainingTop200Response']['response']['rmse'] for k,v in fold['models'].items()}}),flush=True)
    write(a.out/'results.json',results)

if __name__=='__main__':main()
