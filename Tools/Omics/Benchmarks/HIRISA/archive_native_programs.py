#!/usr/bin/env python3
"""Archive completed original-cohort native program scoring and exact replay."""
import argparse,gzip,hashlib,json
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.study/'native-programs'
 complete=read(r/'full-complete.json');check=read(r/'full-independent.json');receipt=read(r/'full/receipt.json');freeze=read(r/'full-output-freeze.json')
 assert complete['status']=='passed' and complete['nativeReplayPassed'] and complete['cells']==1612594
 assert len(complete['commands'])==3 and all(x['returnCode']==0 for x in complete['commands'])
 assert complete['independentCheckSHA256']==sha(r/'full-independent.json') and check['status']=='passed'
 assert check['receiptSHA256']==sha(r/'full/receipt.json') and check['cells']==1612594 and check['programs']==2
 environment=read(r/'name-mapping/runtime-release/environment.json')
 assert complete['binarySHA256']==environment['binarySHA256'] and sha(r/'name-mapping/runtime-release/numivivo')==complete['binarySHA256']
 source=dict(bytes=(r/'full/original.h5ad').stat().st_size,SHA256=sha(r/'full/original.h5ad'))
 assert source['SHA256']==bytes(receipt['source']['bytes']).hex()
 for name,identity in freeze.items():assert identity==dict(bytes=(r/'full'/name).stat().st_size,SHA256=sha(r/'full'/name)),name
 execution=read(r/'full-execution-freeze.json')
 for name,digest in execution['files'].items():assert sha(r/name)==digest,name
 paths=[('bundle/'+name,r/'full'/name) for name in freeze]
 names=['full-complete.json','full-independent.json','full-output-freeze.json','full-execution-freeze.json','full-protocol.json','full-plan.json','run_full.py','check_native_programs.py','check_integration_response.py']
 for phase in ['full-publish','full-independent','full-verify']:names += [phase+'-status.json',phase+'-start.json',phase+'.log']
 paths += [('run/'+name,r/name) for name in names]
 paths += [('runtime/environment.json',r/'name-mapping/runtime-release/environment.json'),('tools/archive_native_programs.py',Path(__file__))]
 reference=a.study/'integration-programs/reference/reference.json';assert sha(reference)==check['referenceSHA256'];paths.append(('reference/reference.json',reference))
 assert len({n for n,_ in paths})==len(paths)
 for _,source_path in paths:assert source_path.is_file() and not source_path.is_symlink()
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for name,source_path in sorted(paths):
  raw=source_path.read_bytes();encoded=gzip.compress(raw,compresslevel=6,mtime=0);target=a.out/(name+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(encoded);digest=hashlib.sha256(raw).hexdigest();assert sha(source_path)==digest
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=digest,storedPath=name+'.gz',storedBytes=len(encoded),storedSHA256=hashlib.sha256(encoded).hexdigest(),gzipEncoded=True))
 manifest=dict(schemaVersion=1,records=records,externalOriginalH5AD=source,externalRuntimeSHA256=complete['binarySHA256'],independentReferenceArchive='2026-09-10-integration-programs',scope='Full native original-cohort program publication, exact replay, every-cell score/detection/total comparison and original metadata identity. All native derived arrays included. Original H5AD and frozen executable remain external with exact identities; independent source/RNA arrays restore from previous integration-programs archive. No biological prediction or preservation gate is promoted by numerical agreement.')
 (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='archived',members=len(records),storedBytes=sum(x['storedBytes'] for x in records),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
