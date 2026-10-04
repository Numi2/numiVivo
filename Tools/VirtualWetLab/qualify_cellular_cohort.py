#!/usr/bin/env python3
"""Executed qualification of admitted sources and bounded observation access."""
import argparse,time
from pathlib import Path
import numpy as np
from wetlab import read,write,require
from measured_cellular_cohort import MeasuredCellularCohort,load,values

def qualify(config,output):
 from prepare_cellular_cohort import row_selection
 reg=read(Path(config).parent/'registration.json')
 started=time.monotonic();owner=MeasuredCellularCohort();a,sources=load(config);checks=[]
 for sid,d in sources.items():
  spec=next(x for x in reg['sources'] if x['id']==sid);metadata=read(Path(config).parent/'source-metadata'/(sid+'.json'));require(row_selection(reg,spec,metadata)==d['rows'],'Deterministic metadata-only cell sampling differs')
  s={'specimen':sid,'population':d['source']['populationID'],'condition':d['source']['conditionID'],'targets':[],'objective':{'genes':['NOC2L'],'preserveGenes':[],'penalty':0}}
  gene=next(g for g,ix in d['features'].items() if len(ix)==1 and not g.startswith('ENSG'));targets=[t for t in d['groups'] if t!='no-intervention'];target=targets[0]
  response=owner.readout(config,None,{'kind':'gene','gene':gene,'target':target,'limit':7,'mode':'measured'},s);arm=response['arms'][0]
  require(len(arm['observedValues'])<=7 and len(arm['controlValues'])<=7,'Readout page bound')
  byid={r['id']:r for r in d['rows']};require(all(byid[i]['observationAccess']=='development-open' for i in arm['observedCellIDs']),'Reserved cells escaped')
  structural_missing=next(g for g in reg['featureAxis'] if g not in d['features'])
  absent=owner.readout(config,None,{'kind':'gene','gene':structural_missing,'target':target,'limit':7},s);require(absent['arms'][0]['observed'] is None and not absent['measured'],'Missing feature became measured zero')
  tile=owner.readout(config,None,{'kind':'expression-tile','genes':[gene,'__UNMEASURED__'],'target':target,'limit':7},s);require(tile['featureMask']==[True,False] and all(v[1] is None for v in tile['values']),'Tile mask differs')
  for q,selection in [({'kind':'gene','gene':gene,'target':target,'mode':'predict'},s),({'kind':'gene','gene':gene,'target':target,'limit':129},s),({'kind':'gene','gene':gene,'target':target},{**s,'condition':'unsupported'})]:
   try:owner.readout(config,None,q,selection)
   except ValueError:pass
   else:raise AssertionError('Unsupported readout accepted')
  reserved=[i for i,r in enumerate(d['rows']) if r['observationAccess']!='development-open']
  if reserved:
   try:values(d,reserved[:1],gene,'raw_counts')
   except ValueError:pass
   else:raise AssertionError('Direct reserved row read accepted')
  plans={task:read(Path(config).parent/sid/(task+'.json')) for task in ('known-target','held-target','held-context')}
  held=set(plans['held-target']['heldOutTargetIDs']);require(not any(x['targetID'] in held for x in plans['held-target']['assignments'] if x['partition']=='training'),'Held target leaked into training')
  held_contexts=set(plans['held-context']['heldOutContextIDs']);require(not any(x['contextID'] in held_contexts for x in plans['held-context']['assignments'] if x['partition']=='training'),'Held context leaked into training')
  require(all(r['guideAssignment'] is None or isinstance(r['guideAssignment'],str) for r in d['rows']),'Guide assignments lost')
  # Metadata partition identities, not resampled bags, determine access.
  require(all(r['control'] or r['observationAccess']!='development-open' or all(v=='training' for v in r['partitions'].values()) for r in d['rows']),'Cross-task access leakage')
  checks.append({'source':sid,'boundedRead':True,'absentFeatureMask':True,'unknownConditionRejected':True,'predictionBlocked':True,'reservedRowsProtected':len(reserved),'guideIdentityRetained':True,'separateTargetContextSplits':True})
 result={'status':'passed','checks':checks,'seconds':time.monotonic()-started,'biologicalPromotion':False,'trainingExecuted':False,'scope':'Data fidelity, partitions, capability and bounded observation access; not biological or interaction qualification'};write(output,result);return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--output',required=True);print(qualify(**vars(p.parse_args())))
