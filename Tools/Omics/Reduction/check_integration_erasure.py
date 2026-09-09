#!/usr/bin/env python3
"""Supplementary control: remove each type-specific condition mean.

The original global-centering control and its failed acceptance gate are retained.
This is a label-informed diagnostic, never an integration candidate or gate retune.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from check_full_integration import metrics
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference',type=Path,required=True);p.add_argument('--protocol',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
prior=json.loads((a.reference/'checks.json').read_text());protocol=json.loads(a.protocol.read_text())
spec=dict(reason='Global condition centering did not meet the original 0.10 accuracy-loss control requirement in mixed-cell-type Kang. Remove cell-type/condition means to represent cell-type-specific responses.',operation='For every source cell-type/condition group, subtract its full-cohort PCA mean. Labels are deliberately used only by this destructive control.',minimumRequiredAccuracyLoss=protocol['gateMargins']['minimumNegativeControlAccuracyLoss'],originalAcceptanceResultsUnchanged=True,originalReferenceSHA256=hashlib.sha256((a.reference/'checks.json').read_bytes()).hexdigest())
(a.out/'design.json').write_text(json.dumps(spec,indent=2)+'\n')
inputs=np.load(a.reference/'evaluation-inputs.npz');x=inputs['scores'];types=inputs['cellTypes'];conditions=inputs['conditions'];erased=x.copy()
for t in np.unique(types):
 for c in np.unique(conditions):
  mask=(types==t)&(conditions==c)
  if mask.any():erased[mask]-=erased[mask].mean(axis=0)
m=metrics(erased,inputs['donors'],conditions,types,inputs['program'],protocol['evaluationNeighbors'],a.out/'neighbors.npz')
loss=prior['baseline']['conditionBalancedAccuracy']-m['conditionBalancedAccuracy']
result=dict(status='supplementary-control-measured',metrics=m,conditionAccuracyLoss=loss,controlDetectsLoss=bool(loss>=spec['minimumRequiredAccuracyLoss']),originalGlobalControlPassed=prior['negativeControlsDetectLoss']['condition'],originalNativePreservationGatesPassed=prior['allNativePreservationGatesPassed'],qualification='Diagnostic improvement only; preserves the failed global control and original native/reference acceptance outcomes. No biological or integration promotion.')
(a.out/'checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='metrics'}),flush=True)
