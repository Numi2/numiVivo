#!/usr/bin/env python3
"""Archive repeated-axis interoperability results; bind large real files externally."""
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
 env=read(r/'environment.json');checks=read(r/'checks-v2/summary.json');unique=read(r/'unique-regression/summary.json');cross=read(r/'cross-host/summary.json')
 assert checks['status']=='passed-H5AD-projection-interoperability' and unique['status']==cross['status']=='passed'
 for c in checks['cases']:
  q=r/'checks-v2'/c['name'];assert sha(q/'original.h5ad')==c['sourceSHA256'] and sha(q/'projected.h5ad')==c['outputSHA256'];assert c['allAnnDataSlotsEqual'] and c['run']['exitCode']==c['verify']['exitCode']==0
 for folder,report in [('unique-regression',unique),('cross-host',cross)]:
  for c in report['cases']:
   assert c['allFourPayloadsExact']
   for name,h in c['files'].items():assert sha(r/folder/c.get('name',c.get('case'))/name)==h
 assert checks['binarySHA256']==unique['binarySHA256']==cross['binarySHA256']==env['binarySHA256']==sha(r/'h5ad-check')
 for n,h in env['compiledSources'].items():assert sha(repo/n)==h,n
 summary=dict(status='qualified-native-repeated-axis-interoperability',positiveCases=len(checks['cases']),rejections=len(checks['negativeCases']),uniquePriorCases=len(unique['cases']),crossHostCases=len(cross['cases']),binarySHA256=env['binarySHA256'],realCases=[c for c in checks['cases'] if c['kind'].startswith(('real','complete-source'))],scope='All supported aligned slots and exact stored values on declared complete prepared real sources; legacy encodings and analytical identity/storage restrictions retained. Not biological or million-cell qualification.')
 write(r/'summary.json',summary)
 selected={};external=[]
 for folder in ('checks','checks-v2','unique-regression','cross-host','identity-projection','identity-control','identity-control-import'):
  for q in (r/folder).rglob('*'):
   if not q.is_file() or '__pycache__' in q.parts:continue
   name=str(q.relative_to(r))
   if q.stat().st_size>4*2**20:external.append(dict(path=name,bytes=q.stat().st_size,SHA256=sha(q)));continue
   selected[name]=q
 for q in r.iterdir():
  if q.is_file() and q.suffix in ('.json','.py','.log','.sha256'):selected[q.name]=q
 for q in owner.iterdir():
  if q.is_file() and q.suffix in ('.py','.md'):selected['tools/'+q.name]=q
 for name in ('Sources/NumiVivoKit/Omics/VivoH5ADProjection.swift','Sources/NumiVivoKit/Omics/VivoH5ADProjectionIO.swift','Tools/Omics/H5AD/Main.swift','Tools/Omics/H5AD/build.sh'):selected['source/'+name]=repo/name
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for name,q in sorted(selected.items()):
  assert not q.is_symlink();raw=q.read_bytes();digest=hashlib.sha256(raw).hexdigest();data=gzip.compress(raw,compresslevel=6,mtime=0);target=a.out/(name+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data);assert sha(q)==digest
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=digest,storedPath=name+'.gz',storedBytes=len(data),storedSHA256=sha(target),gzipEncoded=True))
 write(a.out/'manifest.json',dict(schemaVersion=1,baseCommit=env['baseCommit'],records=records,externalFiles=external,externalRoot=str(r),externalBinary=dict(path=str(r/'h5ad-check'),SHA256=env['binarySHA256']),scope='Complete format fixtures/native artifacts, every real plan/report/receipt and all checks retained. Large original/projected real H5AD and count JSON payloads remain externally retained and hash-bound, with every cell/feature checked. Earlier source-selection and diagnostic-assertion failures preserved.'))
 print(json.dumps(dict(members=len(records),storedBytes=sum(c['storedBytes'] for c in records),externalFiles=len(external),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
