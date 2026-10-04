#!/usr/bin/env python3
"""Post-reveal diagnostics distinguish target information from generic response.

Split cells and capture runs describe sampling/technical variation only.
They are not independent biological accuracy limits.
"""
import argparse
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file
from wetlab import read,write,require
from intervention_design import bound,observations
from prepare_intervention_design import load_source,axis,orthology,normalized,V03

def run(campaign,records,output):
    campaign,inputs,_=bound(campaign);rows=read(inputs/'reserved-rows.json');q=load_file(str(inputs/'reserved.safetensors'));pred=load_file(str(campaign/'all-predictions.safetensors'));features=read(inputs/'features.json');mapping=orthology(V03/'research/MGI-HOM_MouseHumanSequence.rpt');specs=read(inputs/'prepared.json')['sources'];out=[]
    for root in records:
        reg=read(root/'registration.json');comparison=read(root/'comparison.json');p=read(root/'prediction.json');ix=p['rows'];genes=reg['selection']['objective']['genes'];cols=[features.index(g) for g in genes];observation=load_file(str(root/'observations.safetensors'));y=observation['mean'];delta=y-q['context'][ix]
        s=next(s for s in specs if s['id']==rows[ix[0]]['source']);a=load_source(s);ax=axis(a,mapping,s['species']);axiscols=[ax[g] for g in features];rng=np.random.default_rng(271828);split_errors=[];gene_coverage=[];sampling=[]
        control=normalized(a,rows[ix[0]]['controls'],axiscols)
        for j,i in enumerate(ix):
            r=rows[i];values=normalized(a,r['rows'],axiscols);perm=rng.permutation(len(values));first,second=np.array_split(perm,2)
            split_errors.append(float(np.sqrt(np.mean((values[first].mean(0)-values[second].mean(0))**2))))
            gene_coverage.append({'target':r['target'],'objectiveGeneDetectionFraction':(values[:,cols]>0).mean(0).tolist(),'cells':len(values)})
            # Conditional sample-mean SE includes control sampling, never
            # promoted to between-animal or cross-study model uncertainty.
            score=values[:,cols].mean(1);cs=control[:,cols].mean(1)
            sampling.append({'target':r['target'],'conditionalProgramMeanSE':float(np.sqrt(score.var(ddof=1)/len(score)+cs.var(ddof=1)/len(cs)))})
        del a
        models={}
        for model,value in pred.items():
            pd=value[ix]-q['context'][ix];centered=pd-pd.mean(0);truth=delta-delta.mean(0)
            models[model]={'predictedBetweenTargetSD':float(pd[:,cols].mean(1).std()),'observedBetweenTargetSD':float(delta[:,cols].mean(1).std()),
              'centeredResponseRMSE':float(np.sqrt(np.mean((centered-truth)**2))),
              'centroidIdentityAccuracy':float(np.mean(np.argmin(((pd[:,None]-delta[None,:])**2).mean(2),axis=1)==np.arange(len(ix))))}
        distribution={}
        for model in ('target-descriptor','target-ID','shuffled-target-descriptor'):
            native=load_file(str(campaign/model/'reserved-prediction/prediction.safetensors'));v=native['variance'][ix];mu=native['mean'][ix]
            distribution[model]={'marginalGaussianCrossEntropy':float(np.mean(.5*(np.log(2*np.pi*v)+(observation['variance']+(mu-y)**2)/v))),'scope':'Expected marginal Gaussian NLL using measured cell mean and variance; cell dispersion is not transfer uncertainty'}
        out.append({'distributionMetrics':distribution,'context':reg['selection']['context'],'objectiveGenes':genes,'models':models,'medianSplitCellRMSE':float(np.median(split_errors)),
          'objectiveDetection':gene_coverage,'conditionalSampling':sampling,'controlProgramDetection':(control[:,cols]>0).mean(0).tolist(),
          'biologicalReproducibility':'UNAVAILABLE: insufficient independently resolved evaluation units','limits':'Post-reveal explanation of failures; no new promotion metric or threshold. Split-cell agreement is not independent biological accuracy.'})
    write(output,{'contexts':out,'changesPromotionCriteria':False});print(output)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--record',type=Path,action='append',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.campaign,a.record,a.output)
