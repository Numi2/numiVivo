#!/usr/bin/env python3
"""Replay the retained transfer models without rewriting registered paths/seals."""
import argparse,json,hashlib,subprocess
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def replay(campaign,binary,output):
 reg=json.loads((campaign/'registration.json').read_text());assert sha(binary)==reg['binarySHA256'],'Runtime differs from registration';output.mkdir(exist_ok=False);rows=[]
 for fold in reg['folds']:
  for representation in reg['representations']:
   folder=campaign/fold['id']/representation;refit=folder/'refit';seal=json.loads((folder/'prediction-seal.json').read_text());assert sha(folder/'selector.json')==seal['selectorSHA256'],'Selector changed'
   for name,digest in seal['files'].items():
    path=refit/name;assert path.resolve().is_relative_to(refit.resolve()) and sha(path)==digest,'Sealed input changed: '+name
   plan=refit/'plan.json';step=json.loads(plan.read_text())['steps'][0];dest=output/fold['id']/representation;dest.parent.mkdir(exist_ok=True)
   from train_intervention_design import native
   native(binary,'predict',plan,refit/'test-query.safetensors',dest,refit/'training'/f'weights-{step}.safetensors')
   actual=load_file(str(dest/'prediction.safetensors'));expected=load_file(str(refit/'prediction/prediction.safetensors'));assert actual.keys()==expected.keys()
   assert all(actual[k].dtype==expected[k].dtype and actual[k].shape==expected[k].shape and actual[k].tobytes()==expected[k].tobytes() for k in actual),'Named tensor replay differs'
   rows.append({'fold':fold['id'],'representation':representation,'namedTensors':'byte-exact','serializationOrder':'not required to match','biologicalPromotion':False})
 (output/'replay.json').write_text(json.dumps(rows,indent=2));print(f'{len(rows)} native model predictions replayed exactly by named tensor')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();replay(a.campaign,a.binary,a.output)
