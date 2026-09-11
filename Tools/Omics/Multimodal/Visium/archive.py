#!/usr/bin/env python3
"""Retain compact Visium evidence and bind complete external native archives."""
import argparse,gzip,hashlib,json
from pathlib import Path

def digest(b):return hashlib.sha256(b).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);groups={};records=[]
 def add(group,p,label):
  b=p.read_bytes();groups.setdefault(group,[]).append(dict(path=label,bytes=len(b),SHA256=digest(b),rawUTF8=b.decode()))
 for p in sorted(a.study.iterdir()):
  if p.is_file() and p.suffix in ['.json','.log','.txt','.sh','.sha256']:add('execution',p,'study/'+p.name)
 for directory in ['inputs','native','checks','regression']:
  for p in sorted((a.study/directory).rglob('*')):
   if p.is_file() and not p.is_symlink() and p.suffix in ['.json','.log','.csv']:add('execution',p,'study/'+str(p.relative_to(a.study)))
 for p in sorted(Path(__file__).parent.iterdir()):
  if p.is_file():add('recipe',p,'repo/'+str(p.relative_to(a.repo)))
 for name in ['Sources/NumiVivoKit/Omics/VivoMultiAssayVisium.swift','Sources/NumiVivoKit/Omics/VivoMultiAssayIO.swift','Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift','Tools/Omics/H5AD/build.sh','Tests/NumiVivoIntegrationTests/MultiAssayVisiumTests.swift','Tests/NumiVivoIntegrationTests/MultiAssayTests.swift','Tests/NumiVivoIntegrationTests/MultiAssayH5MUImportTests.swift']:
  add('native-source',a.repo/name,'repo/'+name)
 for name,files in sorted(groups.items()):
  raw=(json.dumps(dict(schemaVersion=1,sourceFiles=files),sort_keys=True,separators=(',',':'))+'\n').encode();encoded=gzip.compress(raw,compresslevel=9,mtime=0);stored=name+'.json.gz';(a.out/stored).write_bytes(encoded)
  records.append(dict(sourcePath=name+'.json',sourceBytes=len(raw),sourceSHA256=digest(raw),storedPath=stored,storedBytes=len(encoded),storedSHA256=digest(encoded),gzipEncoded=True,bundledSourceFiles=True))
 fixture=a.study/'regression-inputs/outs/filtered_feature_bc_matrix.h5';b=fixture.read_bytes();encoded=gzip.compress(b,compresslevel=9,mtime=0);(a.out/'regression-counts.h5.gz').write_bytes(encoded)
 records.append(dict(sourcePath='regression-counts.h5',sourceBytes=len(b),sourceSHA256=digest(b),storedPath='regression-counts.h5.gz',storedBytes=len(encoded),storedSHA256=digest(encoded),gzipEncoded=True))
 (a.out/'manifest.json').write_text(json.dumps(dict(schemaVersion=1,records=records),indent=2)+'\n')
 print(json.dumps(dict(groups=len(groups),originalFiles=sum(len(v) for v in groups.values())+1,storedBytes=sum(x['storedBytes'] for x in records))))
if __name__=='__main__':main()
