#!/usr/bin/env python3
"""Retain complete compact transfer evidence and bind full external native bundles."""
import argparse,gzip,hashlib,json
from pathlib import Path

def digest(raw):return hashlib.sha256(raw).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);groups={};records=[]
 def add(group,path,label):
  b=path.read_bytes();groups.setdefault(group,[]).append(dict(path=label,bytes=len(b),SHA256=digest(b),rawUTF8=b.decode()))
 def store(name,b,bundled=False):
  encoded=gzip.compress(b,compresslevel=9,mtime=0);stored=name+'.gz';(a.out/stored).write_bytes(encoded);records.append(dict(sourcePath=name,sourceBytes=len(b),sourceSHA256=digest(b),storedPath=stored,storedBytes=len(encoded),storedSHA256=digest(encoded),gzipEncoded=True,bundledSourceFiles=bundled))
 suffixes={'.json','.log','.sha256','.py','.sh','.swift','.md'}
 for path in sorted(a.study.iterdir()):
  if path.is_file() and path.suffix in suffixes:add('execution',path,'study/'+path.name)
 for name in ['inputs','native','native-qualified','scores','scores-repeat','regression','regression-qualified']:
  for path in sorted((a.study/name).rglob('*')):
   if path.is_file() and not path.is_symlink() and path.suffix in suffixes:add('execution',path,'study/'+str(path.relative_to(a.study)))
 for path in sorted(Path(__file__).resolve().parent.iterdir()):
  if path.is_file():add('recipe',path,'repo/'+str(path.relative_to(a.repo.resolve())))
 for name in ['Sources/NumiVivoKit/Omics/VivoSingleCellReduction.swift','Sources/NumiVivoKit/Omics/VivoSingleCellReference.swift','Sources/NumiVivoKit/Omics/VivoReferenceLogistic.swift','Sources/NumiVivoKit/Omics/VivoSingleCellReferenceIO.swift','Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift','Tests/NumiVivoIntegrationTests/ReferenceFeaturePanelTests.swift','Tests/NumiVivoIntegrationTests/SingleCellReductionTests.swift','Tests/NumiVivoIntegrationTests/SingleCellCohortTests.swift','Tests/NumiVivoIntegrationTests/SingleCellReferenceTests.swift','Tests/NumiVivoIntegrationTests/ReferenceLogisticTests.swift','Tools/Omics/H5AD/build.sh']:
  add('native-source',a.repo/name,'repo/'+name)
 for name,files in sorted(groups.items()):store(name+'.json',(json.dumps(dict(schemaVersion=1,sourceFiles=files),sort_keys=True,separators=(',',':'))+'\n').encode(),True)
 for path in sorted((a.study/'scores').glob('*.npz')):store(path.name,path.read_bytes())
 for path in sorted((a.study/'regression-qualified').glob('*.h5ad')):store('fixture-'+path.name,path.read_bytes())
 (a.out/'manifest.json').write_text(json.dumps(dict(schemaVersion=1,records=records),indent=2)+'\n');print(json.dumps(dict(groups=len(groups),bundledFiles=sum(map(len,groups.values())),storedBytes=sum(r['storedBytes'] for r in records))))
if __name__=='__main__':main()
