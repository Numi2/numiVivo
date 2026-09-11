#!/usr/bin/env python3
"""Retain legacy H5AD qualification and failures, binding large payloads externally."""
import argparse,gzip,hashlib,json
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.root;owner=Path(__file__).parent;repo=owner.parents[3]
 env=read(r/'environment.json');checks=read(r/'checks/summary.json');modern=read(r/'modern-checks/summary.json');native=read(r/'native-checks/summary.json')
 assert checks['status']==native['status']=='passed' and modern['status']=='passed-H5AD-projection-interoperability'
 assert sha(owner/'check_legacy.py')==checks['checkerSHA256'] and sha(owner/'check_legacy_native.py')==native['checkerSHA256']
 for c in checks['cases']:
  q=r/'checks'/c['name'];assert sha(q/'original.h5ad')==c['sourceSHA256'] and sha(q/'projected.h5ad')==c['projectedSHA256']
  assert c['allAnnDataSlotsEqual'] and c['reconstructionPassed']
 for c in native['cases']:
  assert c['allFourPayloadsExact']
  for name,h in c['files'].items():assert sha(r/'checks'/c['name']/name)==h
 for folder in ['unique-regression','repeated-regression']:
  report=read(r/folder/'summary.json');assert report['status']=='passed' and report['binarySHA256']==env['binarySHA256']
  for c in report['cases']:
   assert c['allFourPayloadsExact']
   for name,h in c['files'].items():assert sha(r/folder/c['name']/name)==h
 assert checks['binarySHA256']==modern['binarySHA256']==native['binarySHA256']==env['binarySHA256']==sha(r/'h5ad-check')
 for name,h in env['compiledSources'].items():assert sha(repo/name)==h,name
 cleanup=read(r/'cleanup-plan.json');receipt=read(r/'cleanup-receipt.json');assert receipt['status']=='completed' and receipt['planSHA256']==sha(r/'cleanup-plan.json')
 for record in cleanup['files']:
  retained=Path(record['retainedPath']);assert retained.stat().st_size==record['bytes'] and sha(retained)==record['SHA256']
 summary=dict(status='qualified-native-legacy-H5AD-projection',legacyPositiveCases=len(checks['cases']),legacyControlledRejections=len(checks['rejections']),modernPositiveCases=len(modern['cases']),modernRejections=len(modern['negativeCases']),crossHostCases=len(native['cases']),historicalUniqueCases=10,historicalRepeatedSuiteCases=19,binarySHA256=env['binarySHA256'],realCases=[c for c in checks['cases'] if c['name'].startswith('kang-original')],cleanup=receipt,scope='Complete original Kang counts, annotations and embeddings preserved; exact AnnData semantics for declared legacy representations. No new biological prediction or full product qualification.')
 write(r/'summary.json',summary)
 selected={};external=[]
 for folder in ['checks','checks-v2','checks-v2b','checks-v2c','checks-v5','kang-v5','modern-checks','unique-regression','repeated-regression','native-checks','native-reference']:
  for q in (r/folder).rglob('*'):
   if not q.is_file() or '__pycache__' in q.parts:continue
   name=str(q.relative_to(r))
   if q.stat().st_size>4*2**20:external.append(dict(path=name,bytes=q.stat().st_size,SHA256=sha(q)));continue
   selected[name]=q
 for q in r.iterdir():
  if q.is_file() and q.suffix in ['.json','.py','.log','.sha256','.swift','.txt','.ips']:selected[q.name]=q
 for name in ['check_legacy.py','check_legacy_native.py','archive_legacy.py','check.py','check_unique_replay.py','README.md']:selected['tools/'+name]=owner/name
 for name in ['Sources/NumiVivoKit/Omics/VivoH5ADProjection.swift','Sources/NumiVivoKit/Omics/VivoH5ADLegacyIO.swift','Sources/NumiVivoKit/Omics/VivoH5ADProjectionIO.swift','Tools/Omics/H5AD/build.sh','Tools/Omics/H5AD/Main.swift']:selected['source/'+name]=repo/name
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for name,q in sorted(selected.items()):
  assert not q.is_symlink();raw=q.read_bytes();h=hashlib.sha256(raw).hexdigest();payload=gzip.compress(raw,compresslevel=6,mtime=0);dest=a.out/(name+'.gz');dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(payload);assert sha(q)==h
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=h,storedPath=name+'.gz',storedBytes=len(payload),storedSHA256=sha(dest),gzipEncoded=True))
 binaries=[dict(path=str(q),SHA256=sha(q)) for q in sorted(r.glob('h5ad-check*')) if q.is_file()]
 write(a.out/'manifest.json',dict(schemaVersion=1,baseCommit=env['baseCommit'],records=records,externalRoot=str(r),externalFiles=external,externalExecutables=binaries,scope='Complete format artifacts, all real plans/reports/receipts, exact source identities, reference checks and failed attempts retained. Large complete real payloads are externally hash-bound. Native-host data duplicates were removed only after exact retained-copy and open-handle checks; cleanup-plan.json gives restoration mappings.'))
 print(json.dumps(dict(members=len(records),storedBytes=sum(x['storedBytes'] for x in records),externalFiles=len(external),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
