#!/usr/bin/env python3
"""Archive verified full-cohort Gaussian qualification, retaining pre-repair evidence."""
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
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.root;repo=Path(__file__).resolve().parents[4]
 repair=read(r/'native-phase/repair-comparison.json');assert repair['allNineScoreAndAnchorPayloadsExact'] and repair['allNumericalReportFieldsExact']
 assert sha(r/'cohorts/complete.json')==repair['originalCompleteSHA256'] and sha(r/'native-phase/cohorts/complete.json')==repair['repairedCompleteSHA256']
 for folder in ('cohorts','repaired/cohorts','final/cohorts','qualified/cohorts','release/cohorts','native-phase/cohorts'):
  complete=read(r/folder/'complete.json')
  for name,h in complete['results'].items():
   q=r/folder/name
   if q.exists():assert sha(q)==h,(folder,name)
 for name,h in read(r/'evaluation-complete/complete.json')['files'].items():assert sha(r/'evaluation-complete'/name)==h,name
 assert read(r/'evaluation-complete/execution-freeze.json')['outputsFreezeSHA256']==sha(r/'cohorts/complete.json')
 assert read(r/'native-phase/lifecycle-checks.json')['status']=='passed'
 assert read(r/'kernel-checks-release.json')['status']=='passed' and read(r/'native-phase/kernel-checks.json')['status']=='passed'
 comparison=read(r/'evaluation-complete/original-comparison.json');assert all(not c['metricChanges'] and all(v['changedEntries']==0 for v in c['arrays']) for c in comparison['cohorts'])
 selected={};restoration=[]
 def add(name,p):
  assert name not in selected;selected[name]=p
 for folder in ('cohorts','repaired/cohorts','final/cohorts','qualified/cohorts','release/cohorts','native-phase/cohorts'):
  complete=read(r/folder/'complete.json')
  for name,h in complete['results'].items():
   if name.endswith('.bin'):
    cohort,mode=name.split('/')[0].split('-',1);payload=name.split('/')[-1]
    origin='native-phase/cohorts/'+cohort+'-tiled/scores.bin' if payload=='scores.bin' and mode!='scan' else 'ORIGINAL_COHORT_ROOT/'+cohort+'/mnn/'+payload
    restoration.append(dict(path=folder+'/'+name,SHA256=h,restoreFrom=origin))
  for q in (r/folder).rglob('*'):
   if q.is_file() and q.suffix in ('.json','.log'):add(str(q.relative_to(r)),q)
 for c in ('hagai','kang','ding'):
  q=r/'native-phase/cohorts'/(c+'-tiled')/'scores.bin';assert sha(q)==read(r/'cohorts/complete.json')['results'][c+'-tiled/scores.bin'];add(str(q.relative_to(r)),q)
 for q in (r/'evaluation-complete').rglob('*'):
  if q.is_file():add(str(q.relative_to(r)),q)
 for folder in ('pre-underflow-fix','pre-bridge-extraction','pre-scalar-extraction','pre-native-scalar','pre-native-phase'):
  for q in (r/folder).rglob('*'):
   if q.is_file() and q.suffix in ('.cpp','.h','.swift','.py'):add(str(q.relative_to(r)),q)
 for prefix in ('','repaired/','final/','qualified/','release/','native-phase/'):
  root=r/prefix
  for q in root.iterdir():
   if q.is_file() and q.suffix in ('.json','.log','.py','.swift','.cpp') and q.name not in ('package.json','edit_kernel.py','fix_underflow.py','make_harness.py'):add(str(q.relative_to(r)),q)
  for q in (root/'runtime').glob('*.sha256'):add(str(q.relative_to(r)),q)
  for name in ('plan.json','receipt.json','report.json'):
   q=root/'lifecycle-hagai'/name
   if q.exists():add(str(q.relative_to(r)),q)
 for c in ('hagai','kang','ding'):
  for q in (r/'native-phase'/(c+'-default-threads')).glob('*.json'):add(str(q.relative_to(r)),q)
 for q in (r/'baseline').rglob('*'):
  if q.is_file() and q.suffix in ('.json','.log'):add(str(q.relative_to(r)),q)
 for q in Path(__file__).parent.iterdir():
  if q.is_file():add('tools/'+q.name,q)
 env=read(r/'environment.json')
 for name,h in env['finalSourceSHA256'].items():
  q=repo/name;assert sha(q)==h,name;add('source/'+name,q)
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for name,q in sorted(selected.items()):
  assert q.is_file() and not q.is_symlink() and q.stat().st_size<48*2**20
  raw=q.read_bytes();h=hashlib.sha256(raw).hexdigest();z=gzip.compress(raw,compresslevel=6,mtime=0);target=a.out/(name+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z);assert sha(q)==h
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=h,storedPath=name+'.gz',storedBytes=len(z),storedSHA256=sha(target),gzipEncoded=True))
 write(a.out/'manifest.json',dict(schemaVersion=1,baseCommit=env['baseCommit'],records=records,reusedPayloadRestoration=restoration,externalOriginalInputs=read(r/'cohorts/execution-freeze.json')['inputSHA256'],externalRuntimeBinaries=env['binaries'],scope='Both numerical implementations, failed subnormal check and repair, all complete cohort comparisons, exact replay, complete biological-neighbor/prediction arrays and original-metric comparisons, public lifecycle, tests and source. Stored repaired tiled scores restore identical earlier/replay score payloads; original anchor/scalar payloads require the bound original cohorts. No million-cell or independent biological qualification.'))
 print(json.dumps(dict(members=len(records),storedBytes=sum(x['storedBytes'] for x in records),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
