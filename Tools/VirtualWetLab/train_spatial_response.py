#!/usr/bin/env python3
"""Execute the preregistered native training budget; select only on validation."""
import argparse,subprocess,time,os
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file
from wetlab import read,write,sha,require,timestamp

def execute(root,binary,variants):
 root=Path(root);binary=Path(binary);prepared=read(root/'prepared.json')
 for variant in variants:
  folder=root/variant;plan=read(folder/'plan.json');training=folder/'training'
  receipt={'startedAt':timestamp(),'binarySHA256':sha(binary),'preparedSHA256':sha(root/'prepared.json'),'planSHA256':sha(folder/'plan.json'),'trainSHA256':sha(folder/'chip2.safetensors'),'validationSHA256':sha(folder/'chip3.safetensors'),'preregisteredSteps':plan['steps'],'variant':variant}
  write(folder/'training-invocation.json',receipt)
  start=time.monotonic()
  subprocess.run([str(binary),'spatial-response','train',str(folder/'plan.json'),str(folder/'chip2.safetensors'),str(training)],check=True)
  results=[];validation=load_file(str(folder/'chip3.safetensors'));rows=read(root/'chip3-rows.json')
  for step in plan['steps']:
   out=folder/('validation-'+str(step));weights=training/('weights-'+str(step)+'.safetensors')
   subprocess.run([str(binary),'spatial-response','predict',str(folder/'plan.json'),str(folder/'chip3.safetensors'),str(out),str(weights)],check=True)
   pred=load_file(str(out/'prediction.safetensors'))['mean'];by=[]
   for target in sorted(set(m['target'] for m in rows)):
    idx=[i for i,m in enumerate(rows) if m['target']==target];mask=validation['mask'][idx];err=(pred[idx]-validation['observed'][idx])**2
    by.append({'target':target,'rmse':float(np.sqrt((err*mask).sum()/mask.sum()))})
   results.append({'step':step,'equalTargetRMSE':float(np.mean([x['rmse'] for x in by])),'byTarget':by,'weightsSHA256':sha(weights)})
  chosen=min(results,key=lambda x:(x['equalTargetRMSE'],x['step']))
  write(folder/'selection.json',{'format':'numivivo-spatial-model-selection/v1','selection':'validation only, fixed three-checkpoint budget','candidates':results,'selected':chosen,'seconds':time.monotonic()-start,'biologicalPromotion':False,'testOutcomesRead':False})
  out=folder/'test-prediction';weights=training/('weights-'+str(chosen['step'])+'.safetensors')
  subprocess.run([str(binary),'spatial-response','predict',str(folder/'plan.json'),str(folder/'chip1.safetensors'),str(out),str(weights)],check=True)
  write(folder/'prediction-seal.json',{'format':'numivivo-spatial-prediction-seal/v1','createdAt':timestamp(),'binarySHA256':sha(binary),'weights':str(weights),'weightsSHA256':sha(weights),'inputsSHA256':sha(folder/'chip1.safetensors'),'predictionSHA256':sha(out/'prediction.safetensors'),'selectionSHA256':sha(folder/'selection.json'),'observationPolicy':'held-out outputs not supplied to native prediction','testOutcomesRead':False})
  print(variant,chosen,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--inputs',required=True,type=Path);p.add_argument('--binary',required=True,type=Path);p.add_argument('--variants',nargs='+',default=['no-neighborhood','neighborhood','shuffled-neighborhood','unseen-Cfap410']);a=p.parse_args();execute(a.inputs,a.binary,a.variants)
