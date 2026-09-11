#!/usr/bin/env python3
"""Run unchanged original biological diagnostics on the frozen complete candidate."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from check_full_integration import metrics as donor_metrics,matrix
from evaluate_ding_integration import metrics as ding_metrics

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('root','source','spec','out'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);freeze=read(a.root/'outputs-frozen.json');assert freeze['metricsRead'] is False and freeze['allOriginalCells'] and freeze['allReplayArraysExact'];spec=read(a.spec);owner=Path(__file__).parent.parent
 bindings={}
 for c in ('hagai','kang','ding'):
  inp=Path(spec[c]['inputs']);ref=Path(spec[c]['reference']);old=read(ref/'checks.json');assert sha(inp)==old['inputSHA256'];assert sha(a.root/c/'scores.npz')==freeze['files'][c+'/scores.npz'];bindings[c]=dict(inputsSHA256=sha(inp),referenceChecksSHA256=sha(ref/'checks.json'),scoresSHA256=sha(a.root/c/'scores.npz'))
 margin=read(Path(spec['protocol']))['gateMargins'];write(a.out/'execution-freeze.json',dict(createdUnix=time.time(),outputsFreezeSHA256=sha(a.root/'outputs-frozen.json'),candidateBindings=bindings,protocolSHA256=sha(Path(spec['protocol'])),evaluatorSHA256=sha(Path(__file__)),dependencySHA256={n:sha(owner/n) for n in ('check_full_integration.py','evaluate_ding_integration.py')},margins=margin))
 results=[]
 for c in ('hagai','kang','ding'):
  out=a.out/c;out.mkdir();inputs=np.load(spec[c]['inputs'],allow_pickle=False);old=read(Path(spec[c]['reference'])/'checks.json');baseline=old['baseline'];x=inputs['scores'];np.testing.assert_array_equal(x,matrix(a.source/c/'pca/scores.bin',len(x),x.shape[1]));y=np.load(a.root/c/'scores.npz')['scores'];assert y.shape==x.shape and np.isfinite(y).all();start=time.perf_counter()
  if c=='ding':
   m=ding_metrics(y,inputs,30,4,out,'candidate');g=dict(mixingImproved=m['sameMethodExcess']<baseline['sameMethodExcess'],cellTypesPreserved=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'],everyTypeRecallPreserved=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items()),programsPreserved=all(m['programs'][p]['spearman']>=v['spearman']-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items()),withinStratumProgramsPreserved=all(m['programs'][p]['withinStratumSpearman']>=v['withinStratumSpearman']-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items()))
  else:
   m=donor_metrics(y,inputs['donors'],inputs['conditions'],inputs['cellTypes'],inputs['program'],30,out/'candidate-neighbors.npz');g=dict(mixingImproved=m['sameDonorExcess']<baseline['sameDonorExcess'],conditionPreserved=m['conditionBalancedAccuracy']>=baseline['conditionBalancedAccuracy']-margin['maximumConditionBalancedAccuracyLoss'],programPreserved=m['programSpearman']>=baseline['programSpearman']-margin['maximumProgramSpearmanLoss'],withinStratumProgramPreserved=m['withinStratumProgramSpearman']>=baseline['withinStratumProgramSpearman']-margin['maximumProgramSpearmanLoss'],completeClassifierStrata=not m['missingClassifierStrata'])
   if baseline['cellTypeBalancedAccuracy'] is not None:g.update(cellTypesPreserved=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'],everyTypeRecallPreserved=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items()))
  value=dict(status='complete-original-cohort-measured',cohort=c,cells=len(x),baseline=baseline,metrics=m,gates=g,allMeasuredGatesPassed=all(g.values()),seconds=time.perf_counter()-start,scope='Development approximate matching with full Gaussian correction on inspected full cohorts with unchanged original diagnostics; unavailable/partial labels retained. Not independent biology or prospective prediction.');write(out/'checks.json',value);results.append(dict(cohort=c,cells=len(x),gates=g));print(c,g,flush=True)
 write(a.out/'complete.json',dict(status='completed',results=results,files={str(p.relative_to(a.out)):sha(p) for p in a.out.rglob('*') if p.is_file()}))
if __name__=='__main__':main()
