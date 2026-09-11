#!/usr/bin/env python3
"""Archive the fixed matching intervention, binding all complete step witnesses."""
import argparse,gzip,hashlib,json
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser()
 for n in ('root','out'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();r=a.root;owner=Path(__file__).parent;f=read(r/'candidate/outputs-frozen.json');numerical=read(r/'numerical-checks.json');ev=read(r/'evaluation/complete.json');env=read(r/'environment.json')
 assert f['allOriginalCells'] and f['allReplayArraysExact'] and f['metricsRead'] is False
 assert numerical['outputsFreezeSHA256']==sha(r/'candidate/outputs-frozen.json')
 assert read(r/'evaluation/execution-freeze.json')['outputsFreezeSHA256']==sha(r/'candidate/outputs-frozen.json')
 for n,h in f['files'].items():assert sha(r/'candidate'/n)==h,n
 for n,h in ev['files'].items():assert sha(r/'evaluation'/n)==h,n
 for n,h in read(r/'candidate/freeze.json')['sourceSHA256'].items():assert sha(owner.parent/n)==h,n
 rows=[]
 for c in ('hagai','kang','ding'):
  fit=read(r/'candidate'/c/'report.json');replay=read(r/'candidate'/c/'replay.json');bio=read(r/'evaluation'/c/'checks.json');num=next(v for v in numerical['results'] if v['cohort']==c)
  assert fit['arraySHA256']==replay['arraySHA256']
  rows.append(dict(cohort=c,cells=fit['cells'],fitSeconds=fit['fitSeconds'],replaySeconds=replay['fitSeconds'],anchorPrecision=num['anchorPrecision'],anchorRecall=num['anchorRecall'],anchorGatePassed=num['anchorGatePassed'],assemblyOrderOriginalExact=num['assemblyOrderOriginalExact'],matchingAndAnchorsPreviousExact=num['matchingAndAnchorsPreviousExact'],relativeFrobeniusError=num['relativeFrobeniusError'],maximumKernelDeltaError=num['maximumKernelDeltaError'],gates=bio['gates'],allMeasuredGatesPassed=bio['allMeasuredGatesPassed']))
 write(r/'summary.json',dict(status='completed-fixed-development-intervention',cohorts=rows,scope='Approximate matching held fixed relative to previous local64 trial; full-anchor Gaussian restored. Complete inspected-cohort development evidence, no native API promotion, million-cell or independent biological qualification. Missing/partial labels preserved.'))
 selected={};external=[]
 for folder in ('candidate','evaluation'):
  for q in (r/folder).rglob('*'):
   if not q.is_file() or '__pycache__' in q.parts:continue
   n=str(q.relative_to(r))
   if q.name.startswith('step-') and q.suffix=='.npz':external.append(dict(path=n,bytes=q.stat().st_size,SHA256=sha(q)));continue
   selected[n]=q
 for q in r.iterdir():
  if q.is_file() and q.suffix in ('.json','.log','.py'):selected[q.name]=q
 for q in owner.iterdir():
  if q.is_file():selected['tools/FullGaussianMNN/'+q.name]=q
 for name in ('fit.py','LocalNeighbors.cpp'):
  selected['tools/LocalMNN/'+name]=owner.parent/'LocalMNN'/name
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for n,q in sorted(selected.items()):
  assert not q.is_symlink() and q.stat().st_size<48*2**20
  raw=q.read_bytes();digest=hashlib.sha256(raw).hexdigest();encoded=gzip.compress(raw,compresslevel=6,mtime=0);target=a.out/(n+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(encoded);assert sha(q)==digest
  records.append(dict(sourcePath=n,sourceBytes=len(raw),sourceSHA256=digest,storedPath=n+'.gz',storedBytes=len(encoded),storedSHA256=sha(target),gzipEncoded=True))
 write(a.out/'manifest.json',dict(schemaVersion=1,baseCommit=env['baseCommit'],records=records,externalStepWitnesses=external,externalLibraries=env['libraries'],originalInputs=read(r/'candidate/freeze.json')['inputs'],vendorSHA256=env['vendorSHA256'],scope='All corrected coordinates, matching and mutual anchors, reports, biological neighbor/prediction arrays, independent checks and exact source retained. Complete step snapshots externally retained on both recorded study hosts; no sampling of validation.'))
 print(json.dumps(dict(members=len(records),storedBytes=sum(v['storedBytes'] for v in records),externalWitnessBytes=sum(v['bytes'] for v in external),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
