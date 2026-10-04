#!/usr/bin/env python3
"""Leakage and metric oracles before the grouped development campaign."""
import tempfile
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file
from wetlab import write
from transfer_experiment import preprocessing,data_for,metric

def check():
 with tempfile.TemporaryDirectory() as tmp:
  root=Path(tmp);write(root/'targets.json',{t:{'descriptor':[0.]*149} for t in ['A','B','C']});groups=[{'source':'s','context':'c','target':t,'species':'human','modality':'CRISPRi','unit':'technical'} for t in ['A','B','C']];reg={'groups':groups,'targetManifest':str(root/'targets.json'),'objective':{'genes':['G']}}
  raw={'bags':np.arange(48,dtype=np.float32).reshape(3,8,2),'bagVariance':np.zeros((3,8,2),np.float32),'context':np.ones((3,2),np.float32),'observed':np.asarray([[1,2],[3,4],[5,6]],np.float32)};emb={'A':np.asarray([1,2,3],np.float32),'B':np.asarray([3,2,1],np.float32),'C':np.asarray([99,98,97],np.float32)};var={'s|c':np.asarray([2,1],np.float32)}
  f=root/'a';f.mkdir();a=preprocessing(reg,raw,['G','H'],var,[0,1],'ESM2',emb,f);before=load_file(str(f/'pca.safetensors'));raw['bags'][2]=1e5;raw['observed'][2]=1e5;emb['C'][:]=-1e6;g=root/'b';g.mkdir();b=preprocessing(reg,raw,['G','H'],var,[0,1],'ESM2',emb,g);after=load_file(str(g/'pca.safetensors'));assert all(np.array_equal(before[k],after[k]) for k in before) and np.array_equal(a[-1],b[-1]);assert not np.array_equal(a[1]['C'],b[1]['C'])
  test,rows=data_for(reg,raw,[2],b,'target-ID');assert test['known'].sum()==0
  train,tr=data_for(reg,raw,[0,1],b,'target-ID');assert train['known'].sum()==2
  good=metric(train['observed'],train,tr);bad=metric(2*train['context']-train['observed'],train,tr);assert good['equalGroupRMSE']==0 and bad['equalGroupRMSE']>0
 print('passed: held-target expression/embedding cannot alter fitted PCA or scale; unknown IDs masked; response-error oracle')
if __name__=='__main__':check()
