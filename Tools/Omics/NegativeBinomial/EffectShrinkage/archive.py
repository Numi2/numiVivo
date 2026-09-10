#!/usr/bin/env python3
"""Package checked reports and retained failures; never package a running result."""
import argparse,gzip,hashlib,json,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
assert json.loads((a.root/'check.json').read_text())['status']=='passed-conditional-numerical-checks'
assert json.loads((a.root/'native-complete.json').read_text())['commands']==8
assert '17 tests in 3 suites passed' in (a.root/'focused-tests.log').read_text()
a.out.mkdir(parents=True,exist_ok=False);paths=[]
for name in ['protocol.json','validated-sources.json','native-environment.log','native-runs.json','native-complete.json','native-driver.log','focused-tests.log','reference.log','check.json','restoration.log','restoration-replay.log','restoration-check.json','norman-source-cow.json','remote-report-cow.json','numivivo-report-cow-20260910.json','binary-compression-20260910.json']:
 paths.append(Path(name))
for name in ['kang-default','kang-active','hagai-default','hagai-active']:
 paths.extend(p.relative_to(a.root) for p in (a.root/name).rglob('*') if p.is_file())
for name in ['failure.json','unqualified-reports.json','feature-failures.json.gz','VivoOmicsNegativeBinomial.swift','focused-tests.log','reference.log','count-conversion-reproduction.log','product-build.log','native-runs.json','native-complete.json','native-driver.log','run_native.py']:
 paths.append(Path('attempt-1')/name)
records=[]
for rel in sorted(paths):
 src=a.root/rel;dest=a.out/rel;dest.parent.mkdir(exist_ok=True,parents=True)
 if src.suffix=='.log':
  dest=dest.with_suffix('.log.gz');dest.write_bytes(gzip.compress(src.read_bytes(),mtime=0))
 elif src.suffix=='.gz':subprocess.run(['cp','-c',str(src),str(dest)],check=True)
 else:shutil.copyfile(src,dest)
 raw=dest.read_bytes();assert (gzip.decompress(raw) if src.suffix=='.log' else raw)==src.read_bytes()
 record=dict(path=str(dest.relative_to(a.out)),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))
 if src.suffix=='.log':record.update(logicalPath=str(rel),logicalSHA256=hashlib.sha256(src.read_bytes()).hexdigest())
 records.append(record)
(a.out/'manifest.json').write_text(json.dumps(dict(status='conditional-numerical-qualification-only',entries=records,logicalBytes=sum(r['bytes'] for r in records),retainedFirstAttempt='Full unqualified first-attempt reports remain at paths and hashes in attempt-1/unqualified-reports.json; all original Wald inference was exact. The tested first binary remains on macmini under the run root.',nativeExecutable='macmini:/Users/n/numivivo-effect-shrinkage-20260910/numivivo',sources='External source H5ADs are identified by exact hashes in protocol.json; use restore_bundle.py for archived bundles.'),sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(files=len(records),logicalBytes=sum(r['bytes'] for r in records))))
