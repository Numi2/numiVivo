#!/usr/bin/env python3
"""Development-only matched condition swaps; output movement is not biological validation."""
import argparse
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,timestamp
from train_cohort_recovery import predict

def run(root,out):
    require(not out.exists(),'Retain previous condition probe');out.mkdir();reg=read(root/'registration.json');inputs=Path(reg['inputs']);binary=Path(reg['binary']);a=load_file(str(inputs/'training.safetensors'));rows=read(inputs/'training-rows.json')
    keys=sorted({(r['target'],r['context']) for r in rows if r['source']=='GSE92872'});keys=[k for k in keys if (k[0],'unstimulated' if k[1]=='stimulated' else 'stimulated') in keys]
    idx=[[i for i,r in enumerate(rows) if r['source']=='GSE92872' and (r['target'],r['context'])==k] for k in keys]
    q={'context':np.asarray([a['context'][ix].mean(0) for ix in idx],np.float32),'descriptor':np.asarray([a['descriptor'][ix[0]] for ix in idx]),'target':np.asarray([a['target'][ix[0]] for ix in idx]),'known':np.zeros((len(idx),1),np.float32)}
    actual=np.asarray([(a['observed'][ix]-a['context'][ix]).mean(0) for ix in idx]);swap=np.asarray([keys.index((k[0],'unstimulated' if k[1]=='stimulated' else 'stimulated')) for k in keys]);write(out/'registration.json',{'createdAt':timestamp(),'scope':'exposed GSE92872 development only','keys':keys,'models':['all-mean-source','all-mean-source-RNA-only','all-mean-source-shuffled-metadata'],'checkpoint':5760,'pairedSwap':swap.tolist(),'metrics':'response RMSE against own and paired condition, signed response agreement; paired complete swap must exactly reproduce paired prediction. Mismatched swaps are diagnostics, not real experimental counterfactuals.','binarySHA256':sha(binary),'ownerSHA256':sha(__file__)})
    results={}
    for name in read(out/'registration.json')['models']:
        folder=root/name;weights=folder/'training/weights-5760.safetensors';responses={};queryContext={}
        for mode in ('original','control-swap','descriptor-swap','both-swap'):
            b={k:v.copy() for k,v in q.items()}
            if name.endswith('RNA-only'):b['descriptor'][:,-2:]=0
            if mode in ('control-swap','both-swap'):b['context']=b['context'][swap]
            if mode in ('descriptor-swap','both-swap'):b['descriptor'][:,-2:]=b['descriptor'][swap,-2:]
            path=out/(name+'-'+mode+'.safetensors');save_file(b,str(path));pred,digest=predict(binary,folder/'plan.json',path,weights);responses[mode]=pred['mean']-b['context'];queryContext[mode]=b['context']
            save_file(pred,str(out/(name+'-'+mode+'-prediction.safetensors')))
        results[name]={'weightsSHA256':sha(weights),'pairedCompleteSwapMaxAbs':float(np.max(abs(responses['both-swap']-responses['original'][swap]))),'conditions':[]}
        for i,k in enumerate(keys):
            row={'target':k[0],'condition':k[1],'pairedCondition':keys[swap[i]][1],'observedConditionResponseDistance':float(np.sqrt(np.mean((actual[i]-actual[swap[i]])**2)))}
            for mode,v in responses.items():
                row[mode]={'ownResponseRMSE':float(np.sqrt(np.mean((v[i]-actual[i])**2))),'pairedResponseRMSE':float(np.sqrt(np.mean((v[i]-actual[swap[i]])**2))),'responseMovement':float(np.sqrt(np.mean((v[i]-responses['original'][i])**2)))}
            results[name]['conditions'].append(row)
    write(out/'results.json',results)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.campaign,a.output)
