#!/usr/bin/env python3
"""Archive a rejected exact-tree development trial and bind unchanged large outputs."""
import argparse,gzip,hashlib,json,re
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('root','out'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();r=a.root;c=r/'cohorts-v2';complete=read(c/'complete.json');assert complete['status']=='passed'
 for n,h in complete['results'].items():assert sha(c/n)==h
 env=read(r/'environment.json');assert sha(r/'mnn-tree-check')==env['binaries']['mnn-tree-check']==complete['binarySHA256']
 assert 'passed after' in (r/'tests-v3.log').read_text() and 'with 4 tests' in (r/'tests-v3.log').read_text()
 assert read(r/'tree-checks.json')['status']=='passed' and (r/'sanitizers.log').read_text().startswith('passed:')
 prototype=read(r/'prototype/source-manifest.json')
 for n,h in prototype['files'].items():assert sha(r/'prototype'/n)==h
 inputs=read(c/'execution-freeze.json')['inputSHA256'];restoration={};summary=[]
 for cohort in ('hagai','kang','ding'):
  modes={}
  for mode in ('scan','tree','tree-replay'):
   q=c/(cohort+'-'+mode);checks=read(q/'checks.json');report=read(q/'report.json');assert checks['allOriginalScoreAndAnchorBytesExact'] and checks['allAssemblyStepsExact']
   log=(c/(cohort+'-'+mode+'.log')).read_text();rss=re.search(r'(\d+)\s+maximum resident set size',log)
   modes[mode]=dict(ownerSeconds=checks['seconds'],peakResidentBytes=int(rss[1]),distanceScalarTerms=report['distanceScalarTerms'],kernelScalarTerms=report['kernelScalarTerms'],treeMatching=report.get('treeMatching'))
   for name in ('scores.bin','anchors.bin'):
    f=q/name;h=sha(f);source=cohort+'/mnn/'+name;assert h==inputs[source];restoration[str(f.relative_to(r))]=dict(sourceCohortRelativePath=source,bytes=f.stat().st_size,SHA256=h)
  summary.append(dict(cohort=cohort,cells=read(c/(cohort+'-scan/checks.json'))['cells'],modes=modes))
 write(r/'summary.json',dict(status='completed-not-promoted',cohorts=summary,allOriginalCoordinatesAndAnchorsExact=True,treeMatchingNotPromoted=True,productionChange='Only restore the missing VivoBufferedCountRecords source to the scoped H5AD build list.',scope='Single scan and two tree owner runs per complete original cohort on physical M4 Pro; not repeated end-to-end CLI benchmarking, independent biology or million-cell qualification. Existing biological margins, unavailable strata and partial labels are unchanged because score bytes are exact.'))
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 selected=[]
 for folder in ('cohorts','cohorts-v2','prototype'):
  selected += [(str(f.relative_to(r)),f) for f in (r/folder).rglob('*') if f.is_file() and f.suffix not in ('.bin','.pyc') and '__pycache__' not in f.parts]
 for name in ('environment.json','summary.json','attempts.json','prior-artifacts.json','tree-checks.json','sanitizers.log','SanitizeMain.cpp','TestMain.swift','build-attempt-1.log','build.log','test-build.log','tests.log','tests-v2.log','tests-v3.log','production-restored-build.log'):
  f=r/name
  if f.exists():selected.append((name,f))
 selected += [('runtime/sources.sha256',r/'runtime/sources.sha256'),('production-restored/sources.sha256',r/'production-restored/sources.sha256'),('archive_mnn_tree.py',Path(__file__))]
 for n,f in sorted(selected):
  assert not f.is_symlink() and f.stat().st_size<32*2**20
  raw=f.read_bytes();digest=hashlib.sha256(raw).hexdigest();z=gzip.compress(raw,compresslevel=6,mtime=0);out=a.out/(n+'.gz');out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(z);assert sha(f)==digest
  records.append(dict(sourcePath=n,sourceBytes=len(raw),sourceSHA256=digest,storedPath=n+'.gz',storedBytes=len(z),storedSHA256=sha(out),gzipEncoded=True))
 write(a.out/'manifest.json',dict(schemaVersion=1,records=records,baseCommit=env['baseCommit'],externalExecutables=env['binaries'],sourceBindingsSHA256=sha(r/'environment.json'),externalUnchangedOutputRestoration=restoration,priorCohortManifestSHA256=sha(r/'prior-artifacts.json'),summarySHA256=sha(r/'summary.json'),scope='Rejected development prototype preserved with full modified source overlay. Restore original base commit in an isolated checkout, then overlay prototype files before building. Production API is unchanged. Every omitted score/anchor payload exactly matches its recorded original MNN cohort source; all small results, attempts and checks are retained.'))
 print(dict(members=len(records),storedBytes=sum(v['storedBytes'] for v in records),manifestSHA256=sha(a.out/'manifest.json')))
if __name__=='__main__':main()
