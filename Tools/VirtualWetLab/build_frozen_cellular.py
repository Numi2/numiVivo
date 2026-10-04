#!/usr/bin/env python3
"""Admit an existing frozen transfer checkpoint; no fitting or new source ingestion."""
import argparse,shutil
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,inventory,timestamp

def build(campaign,fold,representation,out):
 campaign=Path(campaign);out=Path(out);out.mkdir(parents=True,exist_ok=False);reg=read(campaign/'registration.json');f=next(x for x in reg['folds'] if x['id']==fold);d=campaign/fold/representation;sel=read(d/'selector.json');seal=read(d/'prediction-seal.json');require(sha(d/'selector.json')==seal['selectorSHA256'],'Selector changed')
 for p,h in seal['files'].items():require(sha(d/'refit'/p)==h,'Frozen checkpoint input changed')
 refit=d/'refit';pre=read(refit/'preprocessing.json');features=pre['features'];allfeatures=read(campaign/'raw/features.json');cols=[allfeatures.index(g) for g in features];raw=load_file(str(campaign/'raw/population.safetensors'));groups=reg['groups'];cases=[]
 for qi,gi in enumerate(f['test']):
  g=groups[gi];cases.append({'index':qi,'groupIndex':gi,'specimenID':g['source'],'populationID':g['source']+'|'+g['context'],'conditionID':g['context'],'target':g['target'],'species':g['species'],'modality':g['modality'],'unit':g['unit'],'biologicalUnitStatus':'source-qualified technical block; independent biological identity unresolved','controlCellIDs':[g['source']+':'+str(i) for i in g['controls']],'outcomeCellIDs':[g['source']+':'+str(i) for i in g['rows']],'guideStatus':'source group assignment; no additional resolved guide identity'})
 baseline=[]
 for gi in f['test']:
  g=groups[gi];ids=[j for j in f['training'] if groups[j]['target']==g['target'] and groups[j]['modality']==g['modality']];scope='target-and-modality matched training response'
  if not ids:ids=[j for j in f['training'] if groups[j]['modality']==g['modality']];scope='training modality mean; target unavailable'
  baseline.append((raw['observed'][ids][:,cols]-raw['context'][ids][:,cols]).mean(0));cases[len(baseline)-1]['baselineScope']=scope
 save_file({'delta':np.asarray(baseline,np.float32)},str(out/'baseline.safetensors'));save_file({'mean':raw['observed'][f['test']][:,cols],'mask':np.ones((len(cases),len(cols)),np.float32)},str(out/'observations.safetensors'))
 for src,name in [(refit/'plan.json','plan.json'),(refit/'test-query.safetensors','query.safetensors'),(refit/'prediction/prediction.safetensors','preview.safetensors'),(refit/'training'/f"weights-{sel['native']['step']}.safetensors",'weights.safetensors'),(d/'selector.json','selector.json'),(d/'prediction-seal.json','source-prediction-seal.json')]:shutil.copy2(src,out/name)
 write(out/'features.json',features);write(out/'cases.json',cases);write(out/'source-registration.json',reg);write(out/'preprocessing.json',pre)
 runtime=out/'runtime';runtime.mkdir();shutil.copy2(reg['binary'],runtime/'numivivo');shutil.copytree(Path(reg['binary']).parent/'mlx-swift_Cmlx.bundle',runtime/'mlx-swift_Cmlx.bundle');require(sha(runtime/'numivivo')==reg['binarySHA256'],'Runtime mismatch')
 specimens=[{'id':s,'title':s,'inferenceSupported':True,'sourceSHA256':next(x['sha256'] for x in reg['sources'] if x['id']==s)} for s in sorted({c['specimenID'] for c in cases})];populations=[{'id':p,'title':p.split('|')[-1],'specimenID':next(c['specimenID'] for c in cases if c['populationID']==p),'conditionID':next(c['conditionID'] for c in cases if c['populationID']==p)} for p in sorted({c['populationID'] for c in cases})]
 artifact={'format':'numivivo-frozen-cellular-assay/v1','family':'learned-cellular-response','id':'cellular-coarse-target0-v04','modelVersion':'v04-transfer-held-target0-coarse-step'+str(sel['native']['step']),'title':'Frozen cellular response · coarse target representation','presentation':'population','geometry':None,'createdAt':timestamp(),'features':features,'featureAxis':'training-fold-selected human symbols; retained historical intersection, not a heterogeneous-panel importer','specimens':specimens,'populations':populations,'conditions':[{'id':c,'title':c} for c in sorted({c['conditionID'] for c in cases})],'targets':[{'target':t} for t in sorted({c['target'] for c in cases})],'runtime':{'binary':'runtime/numivivo','sha256':reg['binarySHA256']},'models':{'native':{'weights':'weights.safetensors','sha256':sha(out/'weights.safetensors')}},'observationAccess':'exposed-development; sealed separately per experiment','biologicalPromotion':False,'evidence':'MODEL INFERENCE','limits':['Frozen one-context-vector cellular model; population-input comparison has not run','Measured population endpoints only; individual cell distributions unavailable in this artifact','Technical development transfer; no independent biological validation','No geometry, receiving-population or preservation support'],'sourceCitation':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90546','scientificOwners':{name:sha(Path(__file__).with_name(name)) for name in ['train_intervention_design.py','train_cohort_recovery.py','design_objective.py']},'files':inventory(out)}
 write(out/'assay.json',artifact);return out/'assay.json'
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--fold',default='held-target-0');p.add_argument('--representation',default='coarse');p.add_argument('--output',type=Path,required=True);a=p.parse_args();print(build(a.campaign,a.fold,a.representation,a.output))
