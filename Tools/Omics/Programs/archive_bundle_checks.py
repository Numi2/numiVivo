#!/usr/bin/env python3
"""Archive native program-bundle qualification, keeping full-cohort execution separate."""
import argparse,gzip,hashlib,json
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.root
 checks=read(r/'release-cli/checks.json');env=read(r/'name-mapping/runtime-release/environment.json')
 assert checks['status']=='passed' and checks['commands']==25 and checks['expectedRejections']==13
 assert checks['binarySHA256']==env['binarySHA256']
 assert read(r/'name-mapping/admission-validation-complete.json')['status']=='passed'
 repo=Path(__file__).resolve().parents[3]
 for name,digest in env['patchFiles'].items():assert sha(repo/name)==digest,name
 paths=[]
 for folder in ['fixtures','debug-cli','release-cli','name-mapping']:
  paths += [('study/'+str(p.relative_to(r)),p) for p in (r/folder).rglob('*') if p.is_file() and p.suffix in {'.json','.log','.py','.stdout','.stderr','.h5ad'}]
 for name in ['initial-admission-failure.json','collection-initial-failure.json','cleanup-result.json','design-decision.json','full-plan.json','full-protocol.json','run_full.py','retire_build_cache.py',
              'admission-tests.log','admission-tests-status.json','admission-release-build.log','admission-release-build-status.json','admission-validation-complete.json','admission-validation-start.json','validate_remote.py']:
  paths.append(('study/'+name,r/name))
 for name in env['patchFiles']:paths.append(('source/'+name,repo/name))
 for name in ['check_bundle.py','archive_bundle_checks.py']:paths.append(('tools/'+name,Path(__file__).with_name(name)))
 cleanup=r.parent/'storage-cleanup/program-build-cache/result.json'
 if cleanup.exists():paths.append(('study/storage-cleanup/result.json',cleanup))
 assert len({name for name,_ in paths})==len(paths)
 for _,source in paths:assert source.is_file() and not source.is_symlink(),source
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for name,source in sorted(paths):
  raw=source.read_bytes();encoded=gzip.compress(raw,compresslevel=6,mtime=0);target=a.out/(name+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(encoded)
  digest=hashlib.sha256(raw).hexdigest();assert sha(source)==digest
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=digest,storedPath=name+'.gz',storedBytes=len(encoded),storedSHA256=hashlib.sha256(encoded).hexdigest(),gzipEncoded=True))
 manifest=dict(schemaVersion=1,records=records,runtimeSHA256=env['binarySHA256'],scope='Twelve native tests in two suites, complete release build and25CLIchecks including13expected rejections; independent controlled counts and unchanged JSON arithmetic. Exact final modified source bytes and full production hashes retained. Initial coordinator root rejection and superseded ID-only build retained. Complete HIRISA protocol is frozen here; its pending native execution is not asserted by this archive. Frozen runtime remains externally stored at native-programs/name-mapping/runtime-release/numivivo.')
 (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='archived',members=len(records),storedBytes=sum(v['storedBytes'] for v in records),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
