#!/usr/bin/env python3
"""Archive complete development results; preserve bulky step witnesses externally."""
import argparse,gzip,hashlib,json
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('root','out'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();r=a.root;owner=Path(__file__).parent;f=read(r/'candidate-v2/outputs-frozen.json');numerical=read(r/'numerical-checks.json');evaluation=read(r/'evaluation/complete.json');env=read(r/'environment.json');assert f['allOriginalCells'] and f['allReplayArraysExact'] and f['metricsRead'] is False
 assert sha(r/'candidate-v2/outputs-frozen.json')==numerical['outputsFreezeSHA256']
 for n,h in f['files'].items():assert sha(r/'candidate-v2'/n)==h
 for n,h in evaluation['files'].items():assert sha(r/'evaluation'/n)==h
 first=read(r/'attempt-1-status.json');assert first['nativeStatus']==2 and first['cohort']=='ding'
 for n,h in first['preservedFileSHA256'].items():assert sha(r/n)==h
 for c in ('hagai','kang'):assert read(r/'candidate'/c/'report.json')['arraySHA256']==read(r/'candidate'/c/'replay.json')['arraySHA256']
 for name,h in read(r/'candidate-v2/freeze.json')['sourceSHA256'].items():assert sha(owner/name)==h
 for name,h in env['libraries'].items():assert sha(r/name)==h
 summary=[]
 for c in ('hagai','kang','ding'):
  fit=read(r/'candidate-v2'/c/'report.json');again=read(r/'candidate-v2'/c/'replay.json');bio=read(r/'evaluation'/c/'checks.json');num=next(x for x in numerical['results'] if x['cohort']==c)
  summary.append(dict(cohort=c,cells=fit['cells'],fitSeconds=fit['fitSeconds'],replaySeconds=again['fitSeconds'],anchorPrecision=num['anchorPrecision'],anchorRecall=num['anchorRecall'],anchorGatePassed=num['anchorGatePassed'],relativeFrobeniusError=num['relativeFrobeniusError'],gates=bio['gates'],allMeasuredGatesPassed=bio['allMeasuredGatesPassed'],maximumKernelDeltaError=num['maximumKernelDeltaError']))
 write(r/'summary.json',dict(status='completed-development-not-promoted',cohorts=summary,firstAttempt='Global filtered matcher failed Ding work admission after complete Hagai/Kang fit and replay; all evidence retained.',scope='One eligible-level design repair with unchanged parameters and no biological scores read before refit. All original cells and fixed evaluation margins retained. Combined approximate matching/local kernel fails preservation in Kang and Ding; no production promotion or million-cell claim.'))
 selected=[];external=[]
 for folder in ('candidate','candidate-v2','evaluation','attempt-1-source'):
  for q in (r/folder).rglob('*'):
   if not q.is_file() or '__pycache__' in q.parts:continue
   n=str(q.relative_to(r))
   if q.suffix=='.dylib' or q.name.startswith('step-') and q.suffix=='.npz':external.append(dict(path=n,bytes=q.stat().st_size,SHA256=sha(q)));continue
   selected.append((n,q))
 for name in ('attempt-1-status.json','environment.json','summary.json','numerical-checks.json','numerical-checks.log','evaluation-spec.json','evaluation.log','alignment-order-comparison.json','sanitizers.log','local-tests.json','native-tests.json','native-tests-v2.json','fit.log','fit-v2.log','native-reference-timings.json','prior-artifacts.json'):
  q=r/name
  if q.exists():selected.append((name,q))
 selected += [('tools/'+q.name,q) for q in owner.iterdir() if q.is_file() and q.suffix in ('.cpp','.py','.md')]
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for n,q in sorted(selected):
  assert not q.is_symlink() and q.stat().st_size<48*2**20
  raw=q.read_bytes();h=hashlib.sha256(raw).hexdigest();z=gzip.compress(raw,compresslevel=6,mtime=0);target=a.out/(n+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z);assert sha(q)==h
  records.append(dict(sourcePath=n,sourceBytes=len(raw),sourceSHA256=h,storedPath=n+'.gz',storedBytes=len(z),storedSHA256=sha(target),gzipEncoded=True))
 write(a.out/'manifest.json',dict(schemaVersion=1,records=records,baseCommit=env['baseCommit'],externalLibraries=env['libraries'],externalStepWitnessesAndOriginalLibrary=external,originalCohortInputs=read(r/'candidate-v2/freeze.json')['inputs'],vendorSHA256=env['vendorSHA256'],summarySHA256=sha(r/'summary.json'),scope='All corrected scores, mutual anchors, matching neighbors, biological-neighbor/prediction payloads, reports, source, both protocols and failed attempts retained. Full per-step source/bias/query/selection/update witnesses remain externally hash-bound at their study-relative paths; deterministic replay verifies their complete arrays. Prototype is not a production API or biological qualification.'))
 print(dict(members=len(records),storedBytes=sum(v['storedBytes'] for v in records),externalWitnessBytes=sum(v['bytes'] for v in external),manifestSHA256=sha(a.out/'manifest.json')))
if __name__=='__main__':main()
