#!/usr/bin/env python3
"""Archive complete original-data analytical qualification and retained failures."""
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
 env=read(r/'environment.json');route=read(r/'route/summary.json');legacy=read(r/'legacy-checks-v2/summary.json');modern=read(r/'modern-regression/summary.json');native=read(r/'native/native-checks.json');ref=read(r/'native-reference.json')
 assert all(x['status']=='passed' for x in [route,legacy,modern,native])
 assert all(x['binarySHA256']==env['binarySHA256'] for x in [route,legacy,modern,native,ref]) and sha(r/'h5ad-check')==env['binarySHA256']
 assert route['checkerSHA256']==sha(owner/'check_legacy_count_route.py') and native['checkerSHA256']==sha(owner/'check_legacy_count_native.py')
 assert native['referenceSHA256']==sha(r/'native-reference.json') and native['allCountStoreAndAggregatePayloadsExact']
 for folder,files in ref['files'].items():
  for name,record in files.items():
   q=r/'route'/folder/name;assert q.stat().st_size==record['bytes'] and sha(q)==record['SHA256']
 for c in legacy['cases']:
  d=r/'legacy-checks-v2'/c['name'];assert c['allAnnDataSlotsEqual'] and c['reconstructionPassed'];assert sha(d/'original.h5ad')==c['sourceSHA256'] and sha(d/'projected.h5ad')==c['projectedSHA256']
 for c in modern['cases']:
  assert c['allFourPayloadsExact']
  for name,h in c['files'].items():assert sha(r/'modern-regression'/c['name']/name)==h
 for name,h in env['compiledSources'].items():assert sha(repo/name)==h,name
 cleanup=read(r/'cleanup-receipt.json');assert cleanup['status']=='completed' and cleanup['planSHA256']==sha(r/'cleanup-plan.json')
 for record in read(r/'cleanup-plan.json')['files']:
  q=Path(record['retainedPath']);assert q.stat().st_size==record['bytes'] and sha(q)==record['SHA256']
 summary=dict(status='qualified-complete-original-Kang-analytical-route',route=route,legacyCases=len(legacy['cases']),legacyRejections=len(legacy['rejections']),modernHistoricalExactCases=len(modern['cases']),native=native,cleanup=cleanup,scope='Original-source projection, annotation, complete counts and pseudobulk only; no new DE/prediction/generalization claim.')
 write(r/'summary.json',summary)
 selected={};external=[]
 for folder in ['route','legacy-checks','legacy-checks-v2','modern-regression','native']:
  for q in (r/folder).rglob('*'):
   if not q.is_file() or '__pycache__' in q.parts:continue
   name=str(q.relative_to(r))
   if q.suffix in ['.h5ad','.bin'] and q.stat().st_size>4*2**20:external.append(dict(path=name,bytes=q.stat().st_size,SHA256=sha(q)));continue
   selected[name]=q
 for q in r.iterdir():
  if q.is_file() and q.suffix in ['.json','.py','.log','.stdout','.sha256']:selected[q.name]=q
 for name in ['check_legacy_count_route.py','check_legacy_count_native.py','check_legacy.py','check_unique_replay.py','archive_legacy_count_route.py','LEGACY_COUNT_ROUTE.md']:selected['tools/'+name]=owner/name
 for name in ['Sources/NumiVivoKit/Omics/VivoH5ADProjection.swift','Sources/NumiVivoKit/Omics/VivoH5ADCountStore.swift','Sources/NumiVivoKit/Omics/VivoH5ADCountAccess.swift','Sources/NumiVivoKit/Omics/VivoH5ADPseudobulk.swift','Sources/NumiVivoKit/Omics/VivoH5ADAnnotations.swift','Tools/Omics/H5AD/Main.swift','Tools/Omics/H5AD/build.sh']:selected['source/'+name]=repo/name
 failed=read(r/'annotation-before-status.json');assert failed['returncode']==65
 dependencies=[]
 for path in [failed['command'][0],failed['command'][2],'/Users/home/numivivo-file-integration-20260909/prepared/kang-fit.json']:
  q=Path(path);dependencies.append(dict(path=path,bytes=q.stat().st_size,SHA256=sha(q)))
 selected['reference/frozen-kang-fit.json']=Path(dependencies[-1]['path'])
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for name,q in sorted(selected.items()):
  assert not q.is_symlink();raw=q.read_bytes();h=hashlib.sha256(raw).hexdigest();data=gzip.compress(raw,compresslevel=6,mtime=0);dest=a.out/(name+'.gz');dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data);assert sha(q)==h
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=h,storedPath=name+'.gz',storedBytes=len(data),storedSHA256=sha(dest),gzipEncoded=True))
 write(a.out/'manifest.json',dict(schemaVersion=1,baseCommit=env['baseCommit'],records=records,externalRoot=str(r),externalFiles=external,externalDependencies=dependencies,externalBinary=dict(path=str(r/'h5ad-check'),SHA256=sha(r/'h5ad-check')),scope='All complete count/aggregate reports, metadata, plans, receipts, source and failures retained. Large H5AD/count payloads are externally retained and hash-bound. Native duplicate cleanup preserves exact restoration mappings.'))
 print(json.dumps(dict(members=len(records),storedBytes=sum(x['storedBytes'] for x in records),externalFiles=len(external),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
