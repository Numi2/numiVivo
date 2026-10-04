#!/usr/bin/env python3
"""Execute native cellular adapter qualification on an admitted frozen artifact."""
import argparse,tempfile,json
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file
import laboratory
from wetlab import write,sha,read,require

def qualify(config,out):
 a=read(config);owner=laboratory.adapter_for_config(config);cases=read(config.parent/'cases.json');c=next(c for c in cases if c['specimenID']=='GSE90546');targets=[x['target'] for x in cases if x['populationID']==c['populationID']][:2]
 selection={'specimen':c['specimenID'],'population':c['populationID'],'condition':c['conditionID'],'targets':targets,'objective':{'genes':['FOS','JUN'],'preserveGenes':[],'penalty':0}}
 rejected=[]
 for bad in [{**selection,'condition':'invented-condition'},{**selection,'targets':['invented-target']},{**selection,'objective':{'genes':['UNMEASURED'],'preserveGenes':[],'penalty':0}}]:
  try:owner.validate_selection(config,bad)
  except ValueError:rejected.append(True)
  else:raise AssertionError('Unsupported input accepted')
 with tempfile.TemporaryDirectory(prefix='qualify-cellular-') as d:
  run=owner.predict(config,{},Path(d),selection);reg=read(run/'registration.json');p=load_file(str(run/'prediction/prediction.safetensors'));preview=load_file(str(config.parent/'preview.safetensors'));indices=[c['index'] for c in reg['cases']]
  require(all(np.array_equal(p[k],preview[k][indices]) for k in p),'Selected native prediction differs from frozen artifact')
  before=owner.readout(config,run,{'gene':'FOS'},selection);require(all(x['observed'] is None for x in before['arms']),'Inspection exposed observations')
  owner.verify(run,{})
  # Explicitly authorized qualification reveal of already-exposed development data.
  owner.reveal(run,{});evaluation=owner.evaluate_objective(config,run,selection['objective']);verification=owner.verify(run,{})
  relocated=Path(d)/'relocated';run.rename(relocated);require(owner.verify(relocated,{})['exactPredictionReplay'],'Relocated replay failed')
  seal=sha(relocated/'seal.json');(relocated/'query.safetensors').write_bytes(b'corrupt')
  try:owner.verify(relocated,{})
  except ValueError:tamper=True
  else:raise AssertionError('Tampered prediction accepted')
 result={'artifactSHA256':sha(config),'unsupportedInputsRejected':len(rejected),'selectedQueryMatchesFrozenPrediction':True,'inspectionDoesNotReveal':True,'authorizedRevealAndEvaluation':True,'exactNativeReplay':verification,'relocationWithoutOriginal':True,'tamperRejected':tamper,'objectiveEvaluation':evaluation,'softwareQualification':True,'biologicalPromotion':False}
 write(out,result);print(json.dumps({k:v for k,v in result.items() if k!='objectiveEvaluation'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--assay',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();qualify(a.assay,a.output)
