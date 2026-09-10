#!/usr/bin/env python3
"""Archive post-result decoder calibration without replacing frozen failures."""
import argparse,gzip,hashlib,json
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.root
 complete=read(r/'results/complete.json');oracles=read(r/'dense-oracles-complete.json');protocol=read(r/'protocol.json');freeze=read(r/'execution-freeze.json')
 assert complete['status']=='completed-development-evaluation' and len(complete['results'])==7
 assert oracles['status']=='passed' and len(oracles['commands'])==5 and all(v['returnCode']==0 for v in oracles['commands'])
 for name,digest in freeze['files'].items():assert sha(r/name)==digest,name
 for name,digest in complete['results'].items():assert sha(r/'results'/(name+'.json'))==digest
 for name in protocol['representations']:
  check=read(r/(name+'-dense-oracle.json'));assert check['status']=='passed' and check['folds']==79 and check['programs']==2
  assert check['resultSHA256']==complete['results'][name]
 previous=read(Path(__file__).parent/'evidence/2026-09-10-integration-programs/manifest.json')
 previous_map={v['sourceSHA256']:v['sourcePath'] for v in previous['records']}
 restoration={name:dict(archive='2026-09-10-integration-programs',sourcePath=previous_map[digest],SHA256=digest) for name,digest in protocol['previousResultsSHA256'].items()}
 paths=[('study/'+str(p.relative_to(r)),p) for p in r.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix in {'.py','.json','.log'}]
 paths.append(('tools/archive_program_calibration.py',Path(__file__)))
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for name,source in sorted(paths):
  assert not source.is_symlink();raw=source.read_bytes();encoded=gzip.compress(raw,compresslevel=6,mtime=0);target=a.out/(name+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(encoded);digest=hashlib.sha256(raw).hexdigest();assert sha(source)==digest
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=digest,storedPath=name+'.gz',storedBytes=len(encoded),storedSHA256=hashlib.sha256(encoded).hexdigest(),gzipEncoded=True))
 manifest=dict(schemaVersion=1,records=records,previousResultRestoration=restoration,scope='Post-result within-library fitting development, seven original-cohort diagnostics and five independent every-fold centered per-cell SVD checks. Reuses exact prior full-source moments and matrices; preserves original mixed-objective failures.28of32controls sensitive,4insufficient; no complete biological preservation or prospective prediction qualification.')
 (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='archived',members=len(records),storedBytes=sum(x['storedBytes'] for x in records),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
