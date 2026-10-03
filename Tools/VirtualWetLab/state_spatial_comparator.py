#!/usr/bin/env python3
"""Run pinned upstream State PerturbMean on matched spatial training inputs.

External research comparator only; no upstream code or pretrained State weights
are incorporated into the native product. State's CC-BY-NC-SA-4.0 terms apply to
its downloaded source. Inputs are mouse log1p(CPM) anchor-distribution means,
not State's published human benchmark. No hardware-dependent CUDA path is used.
"""
import argparse,sys,urllib.request,hashlib
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from wetlab import write,read,sha,timestamp,require
REV='9bbfe78a434a55205e4de834e1ea99f85f7a3add'
def run(root,research):
 root=Path(root);research=Path(research);pkg=research/'state_comparator';pkg.mkdir(exist_ok=True);(pkg/'__init__.py').touch()
 for name in ('base.py','utils.py','perturb_mean.py'):
  dest=pkg/name
  if not dest.exists():urllib.request.urlretrieve('https://raw.githubusercontent.com/ArcInstitute/state/'+REV+'/src/state/tx/models/'+name,dest)
 out=root/'external-State-PerturbMean';out.mkdir(exist_ok=False)
 write(out/'admission.json',{'revision':REV,'method':'State PerturbMean','license':'CC-BY-NC-SA-4.0; retained outside native product for research comparison','species':'species-agnostic refit on mouse gene-symbol axis; no pretrained human State weights','assay':'same source-bound endpoint log1p(CPM) distributions as native learner','hardware':'Apple CPU; torch 2.8, no CUDA','source':{p.name:sha(p) for p in pkg.glob('*.py')},'limitations':'baseline comparator; not a reproduction of full State ST benchmark scores'})
 sys.path.insert(0,str(research));import torch
 from state_comparator.perturb_mean import PerturbMeanPerturbationModel
 from types import SimpleNamespace
 arrays=load_file(str(root/'neighborhood/chip2.safetensors'));meta=read(root/'chip2-rows.json');g=arrays['context'].shape[1]
 model=PerturbMeanPerturbationModel(input_dim=g,hidden_dim=64,output_dim=g,pert_dim=len(read(root/'prepared.json')['targets']),control_pert='mSafe',gene_decoder_bool=False)
 # Preserve direct/neighbor matching; upstream averages each perturbation's
 # offsets across the supplied cell-type strata. Identical cells are not
 # promoted into biological replicates.
 batches=[]
 for start in range(0,len(meta),64):
  rows=meta[start:start+64];n=len(rows);y=np.r_[arrays['observed'][start:start+n],arrays['context'][start:start+n]]
  types=[r['role']+'|'+r['cellType'] for r in rows]
  batches.append({'pert_cell_emb':torch.from_numpy(y),'pert_name':[r['target'] for r in rows]+['mSafe']*n,'cell_type':types+types})
 model.trainer=SimpleNamespace(datamodule=SimpleNamespace(train_dataloader=lambda:batches));model.on_fit_start();model.eval()
 query=load_file(str(root/'neighborhood/chip1.safetensors'));test=read(root/'chip1-rows.json')
 with torch.no_grad():pred=model({'ctrl_cell_emb':torch.from_numpy(query['context']),'pert_name':[r['target'] for r in test]}).numpy()
 require(np.isfinite(pred).all(),'Nonfinite external prediction');save_file({'mean':pred},str(out/'prediction.safetensors'))
 write(out/'prediction-seal.json',{'createdAt':timestamp(),'predictionSHA256':sha(out/'prediction.safetensors'),'admissionSHA256':sha(out/'admission.json'),'trainSHA256':sha(root/'neighborhood/chip2.safetensors'),'querySHA256':sha(root/'neighborhood/chip1.safetensors'),'testOutcomesRead':False})
 print('Pinned State PerturbMean executed and sealed',pred.shape,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--research',type=Path,required=True);a=p.parse_args();run(a.inputs,a.research)
