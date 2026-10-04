#!/usr/bin/env python3
"""Qualify engineering fit separately from exposed and reserved biological results."""
import argparse
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file
from wetlab import read,write,require,sha
from train_intervention_design import native

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root;d=root/'diagnosis-v2';result=read(d/'results.json');x=load_file(str(d/'training.safetensors'));checks={}
require(len(set(x['target'].tolist()))==3,'Targets collapsed');require(len(np.unique(x['descriptor'][:,:149],axis=0))==3,'Descriptors collapsed');require(np.all(x['known']==0),'Descriptor comparison accidentally used ID');require(np.all(x['mask']==1),'Unexpected missing diagnostic features')
for arm,r in result['results'].items():
 if arm=='descriptor-regression':continue
 rows=read(d/arm/'training/training-diagnostics.json');require(all(q['maskedObservationInvariant'] for q in rows),'Mask leakage')
 require(all(q['measuredLossElements']==16*515 for q in rows),'Loss normalization denominator changed')
 err=max(abs(q['gradientCheck']['analytic']-q['gradientCheck']['finiteDifference']) for q in rows)
 require(err<.0005,'Float32 finite-difference check failed')
 require(all(q['parameterUpdateNorms']['targetResponsePrior.weight']==0 for q in rows),'Frozen prior updated')
 require(rows[2]['parameterUpdateNorms']['descriptorProjection.weight']>0,'Encoder did not update')
 if arm=='mean-adam-ID':require(rows[2]['gradientNorms']['targetEmbedding.weight']>0,'Target ID did not affect gradient')
 checks[arm]={'finiteDifferenceMaxAbsError':err,'maskedOutcomesInvariant':True,'frozenPriorUnchanged':True,'fitGate':r['fitGate']}
require(not checks['mean-sgd']['fitGate'] and checks['mean-adam']['fitGate'],'Optimizer contrast not demonstrated')
require(checks['distribution-adam']['fitGate'] and checks['mean-adam-regularized']['fitGate'],'Restored objective or regularization failed')
# A new process must reproduce the retained native distribution tensors.
arm=d/'distribution-adam';dest=root/'diagnostic-replay';native(root/'runtime/numivivo','predict',arm/'plan.json',arm/'training.safetensors',dest,arm/'training/weights-1440.safetensors');original=load_file(str(arm/'fit-1440/prediction.safetensors'));replay=load_file(str(dest/'prediction.safetensors'));require(all(np.array_equal(v,replay[k]) for k,v in original.items()),'Native replay mismatch')
write(root/'training-qualification.json',{'checks':checks,'exactReplay':True,'architectureChanged':False,'engineeringFitQualified':True,'biologicalPromotion':False,'diagnosisSHA256':sha(d/'results.json'),'finding':'AdamW recovers target-specific fit at fixed architecture, data and budget. Removing variance or decay alone does not; response scaling is unnecessary. Does not identify a unique optimal optimizer or establish biological transfer.'})
