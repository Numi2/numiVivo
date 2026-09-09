#!/usr/bin/env python3
"""Evaluate declared adaptive candidates with unchanged frozen-cohort margins."""
import argparse,hashlib,json,platform,time
from importlib.metadata import version
from pathlib import Path
import harmonypy
import numpy as np
import pandas as pd
from check_full_integration import independent,matrix,metrics
p=argparse.ArgumentParser(description=__doc__)
for name in ['native','frozen','protocol','margins','out']:p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--cohort',choices=['kang','hagai'],action='append')
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
protocol=json.loads(a.protocol.read_text());original=json.loads(a.margins.read_text());margin=original['gateMargins']
assert version('harmonypy')=='2.0.0'
for cohort in a.cohort or ['kang','hagai']:
 out=a.out/cohort;out.mkdir();prior=json.loads((a.frozen/f'reference-{cohort}/checks.json').read_text())
 inputs=np.load(a.frozen/f'reference-{cohort}/evaluation-inputs.npz');x=inputs['scores'];donors=inputs['donors'];conditions=inputs['conditions'];types=inputs['cellTypes'];program=inputs['program']
 np.testing.assert_array_equal(x,matrix(a.native/cohort/'pca/scores.bin',*x.shape))
 metadata=json.loads((a.native/cohort/'pca/metadata.json').read_text());oldmeta=json.loads((a.frozen/f'real-evidence/{cohort}/pca/metadata.json').read_text());assert metadata==oldmeta
 baseline=prior['baseline']
 result=dict(status='in-progress',cohort=cohort,cells=len(x),protocolSHA256=hashlib.sha256(a.protocol.read_bytes()).hexdigest(),frozenChecksSHA256=hashlib.sha256((a.frozen/f'reference-{cohort}/checks.json').read_bytes()).hexdigest(),versions={k:version(k) for k in ['numpy','scipy','scikit-learn','harmonypy']},platform=platform.platform(),baseline=baseline,originalNegativeControls=prior['negativeControlsDetectLoss'],native=[],references=[])
 def save():(out/'checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 def measure(label,scores):
  start=time.perf_counter();m=metrics(scores,donors,conditions,types,program,original['evaluationNeighbors'],out/(label+'-neighbors.npz'));m['evaluationSeconds']=time.perf_counter()-start;print(cohort+' '+label+' measured',flush=True);return m
 def gates(m):
  g=dict(mixingImproved=m['sameDonorExcess']<baseline['sameDonorExcess'],conditionPreserved=m['conditionBalancedAccuracy']>=baseline['conditionBalancedAccuracy']-margin['maximumConditionBalancedAccuracyLoss'],programPreserved=m['programSpearman']>=baseline['programSpearman']-margin['maximumProgramSpearmanLoss'],withinStratumProgramPreserved=m['withinStratumProgramSpearman']>=baseline['withinStratumProgramSpearman']-margin['maximumProgramSpearmanLoss'],completeClassifierStrata=not m['missingClassifierStrata'])
  if baseline['cellTypeBalancedAccuracy'] is not None:
   g['cellTypesPreserved']=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'];g['everyTypeRecallPreserved']=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items())
  return g
 save()
 for seed in protocol['seeds']:
  y,math=independent(a.native/cohort/f'seed-{seed}',x,donors);m=measure(f'native-{seed}',y);result['native'].append(dict(seed=seed,metrics=m,gates=gates(m),independentReconstruction=math));save()
 opts=dict(theta=2,lamb=None,sigma=.1,nclust=100,tau=0,block_size=.05,max_iter_harmony=10,max_iter_kmeans=4,epsilon_cluster=.001,epsilon_harmony=.01,alpha=protocol['alpha'],batch_prop_cutoff=1e-5,ncores=1)
 for seed in protocol['seeds']:
  start=time.perf_counter();h=harmonypy.run_harmony(x,pd.DataFrame({'donor':donors}),['donor'],**opts,random_state=seed,verbose=False);seconds=time.perf_counter()-start;y=np.asarray(h.Z_corr);assert y.shape==x.shape and np.isfinite(y).all();np.savez_compressed(out/f'harmony-{seed}.npz',scores=y)
  m=measure(f'harmony-{seed}',y);result['references'].append(dict(seed=seed,metrics=m,gates=gates(m),integrationSeconds=seconds,objectives=list(h.objective_harmony)));save()
 result['status']='full-cohort-adaptive-native-and-reference-measured';result['allNativePreservationGatesPassed']=all(all(r['gates'].values()) for r in result['native'])
 result['qualification']='Retrospective cohorts, unchanged margins; fixed failures retained. Baseline metrics reused only after all PCA bits and metadata match. Independent source-label validation is pending. No prospective prediction or cross-host speed comparison.';save()
