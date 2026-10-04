#!/usr/bin/env python3
"""Native plan admission and historical inference regression."""
import argparse,json,subprocess,tempfile
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file
from wetlab import read,write,sha,require

def run(binary,diagnosis,out):
    require(not out.exists(),'Retain qualification');out.mkdir();base=read(diagnosis/'mean-adam/plan.json');cases=[('legacy',[240,720,1440],None,True),('unregistered',[1,2],None,False),('registered',[1,2],2,True),('duplicate',[1,1],1,False),('descending',[2,1],2,False),('zero',[0,2],2,False),('over-budget',[1,3],2,False),('unused-budget',[1,2],3,False),('unbounded',[100001],100001,False),('empty',[],2,False),('too-many',list(range(1,34)),33,False)]
    result=[]
    for name,steps,budget,expected in cases:
        p={**base,'steps':steps}
        if budget is not None:p['trainingBudget']=budget
        path=out/(name+'.json');write(path,p)
        r=subprocess.run([str(binary),'spatial-response','train',str(path),str(out/'deliberately-missing-input.safetensors'),str(out/(name+'-output'))],capture_output=True,text=True)
        rejected='invalid spatial learning plan' in r.stdout+r.stderr
        require(rejected!=expected,'Unexpected admission '+name);result.append({'name':name,'admittedBeforeInputRead':not rejected,'expected':expected,'exit':r.returncode})
    from train_cohort_recovery import predict
    old=diagnosis/'mean-adam';original=load_file(str(old/'training-prediction/prediction.safetensors')) if (old/'training-prediction/prediction.safetensors').exists() else None
    write(out/'schedule-results.json',{'binarySHA256':sha(binary),'cases':result,'legacyNumerics':'separate exact replay required; admission alone is not replay'})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--diagnosis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.binary,a.diagnosis,a.output)
