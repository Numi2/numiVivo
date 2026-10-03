#!/usr/bin/env python3
"""Full real-source scientific lifecycle and tamper qualification."""
import argparse,json
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file
from learned_spatial import LearnedSpatialResponseAdapter
from wetlab import write,read,require,sha

def check(config,workspace,output):
 adapter=LearnedSpatialResponseAdapter();selection={'specimen':'chip1','population':'neuron|0','targets':['Cfap410','Fasn','Gfap']};checks=[]
 try:adapter.compile(config,{**selection,'targets':['unsupported-target']});raise AssertionError('unsupported accepted')
 except ValueError:checks.append('unsupported intervention rejected')
 run=adapter.predict(config,{},workspace,selection);before=adapter.feature(config,'chip1','Gfap',run);require(all(x['observed'] is None and x['observedRows']==[] for x in before['arms']),'Observation leak before reveal');checks.append('observations unavailable before reveal')
 a=read(config);mask=np.load(config.parent/'model-feature-mask.npy');features=read(config.parent/'features.json');missing=next((features[i] for i,v in enumerate(mask) if not v),None)
 if missing:
  field=adapter.feature(config,'chip1',missing,run);require(all(x['predicted'] is None and x['variance'] is None and x['evidence']=='UNAVAILABLE' for x in field['arms']),'Untrained marker presented as a prediction');checks.append('unmodeled marker remains unavailable')
 require(len(before['arms'])==6 and len(set(x['control'] for x in before['arms']))==1,'Arms do not share reference');checks.append('three interventions and two response roles share one control reference')
 query=run/'neighborhood-query.safetensors';original=query.read_bytes();bad=bytearray(original);bad[-1]^=1;query.write_bytes(bad)
 try:adapter.summary(run);raise AssertionError('tampered query accepted')
 except ValueError:checks.append('changed sealed query rejected')
 finally:query.write_bytes(original)
 adapter.reveal(run,{});after=adapter.feature(config,'chip1','Gfap',run);require(any(x['observed'] is not None for x in after['arms']),'No observations revealed');require(all(x['residual'] is None or abs(x['residual']-(x['predicted']-x['observed']))<1e-6 for x in after['arms']),'Residual sign wrong');checks.append('measured outcomes and correctly signed residuals revealed')
 comparison=run/'comparison.json';original=comparison.read_bytes();changed=read(comparison);changed['biologicalPromotion']=True;comparison.write_text(json.dumps(changed))
 try:adapter.summary(run);raise AssertionError('altered biological verdict accepted')
 except ValueError:checks.append('changed evaluation verdict rejected')
 finally:comparison.write_bytes(original)
 replay=adapter.verify(run,{});checks.append('native mean and variance bit-exact replay; evaluation reconstructs')
 require(read(comparison)['biologicalPromotion'] is False,'Software promoted biology');checks.append('insufficient controls block biological promotion')
 write(output,{'status':'passed','run':str(run),'configSHA256':sha(config),'checks':checks,'replay':replay,'softwareReleaseQualified':True,'biologicalPromotion':False});print(run)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--workspace',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();check(a.config,a.workspace,a.output)
