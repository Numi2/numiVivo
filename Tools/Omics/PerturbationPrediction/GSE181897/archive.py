#!/usr/bin/env python3
"""Archive compact source admission and executed checks; bind full external artifacts."""
import argparse,gzip,hashlib,json
from pathlib import Path

def digest(b):return hashlib.sha256(b).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);groups={}
 def add(group,path,label):
  b=path.read_bytes();groups.setdefault(group,[]).append(dict(path=label,bytes=len(b),SHA256=digest(b),rawUTF8=b.decode()))
 # Retain our audit records and failures. Original author documents and large
 # count/metadata arrays remain at explicitly hashed external locations.
 for p in sorted(a.study.iterdir()):
  if p.is_file() and p.suffix in ['.json','.log'] and p.name!='metadata-values.log':add('admission-and-checks',p,'study/'+p.name)
 for p in sorted((a.study/'handoff').glob('*.json')):add('frozen-handoff',p,'study/handoff/'+p.name)
 for directory in ['native-handoff','native-handoff-qualified']:
  for p in sorted((a.study/directory).rglob('*')):
   if p.is_file() and p.suffix in ['.json','.log'] and p.name!='report.json':add('native-execution',p,'study/'+str(p.relative_to(a.study)))
 for p in sorted(a.study.glob('*.py')):
  if p.name!='run_native.py':add('executed-scripts',p,'study/'+p.name)
 for p in sorted(Path(__file__).parent.iterdir()):
  if p.is_file():add('published-recipe',p,'repo/'+str(p.relative_to(a.repo)))
 records=[]
 for name,files in sorted(groups.items()):
  raw=(json.dumps(dict(schemaVersion=1,sourceFiles=files),sort_keys=True,separators=(',',':'))+'\n').encode();encoded=gzip.compress(raw,compresslevel=9,mtime=0);stored=name+'.json.gz';(a.out/stored).write_bytes(encoded)
  records.append(dict(sourcePath=name+'.json',sourceBytes=len(raw),sourceSHA256=digest(raw),storedPath=stored,storedBytes=len(encoded),storedSHA256=digest(encoded),gzipEncoded=True,bundledSourceFiles=True))
 (a.out/'manifest.json').write_text(json.dumps(dict(schemaVersion=1,records=records),indent=2)+'\n')
 print(json.dumps(dict(groups=len(groups),originalFiles=sum(len(v) for v in groups.values()),storedBytes=sum(x['storedBytes'] for x in records))))
if __name__=='__main__':main()
