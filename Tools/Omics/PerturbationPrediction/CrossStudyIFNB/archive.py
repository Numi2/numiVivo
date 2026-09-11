#!/usr/bin/env python3
"""Archive compact exact-byte evidence; keep large payload hashes explicit."""
import argparse,gzip,hashlib,json
from pathlib import Path

def digest(b):return hashlib.sha256(b).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 groups={}
 def add(group,path,label):
  b=path.read_bytes();groups.setdefault(group,[]).append(dict(path=label,bytes=len(b),SHA256=digest(b),rawUTF8=b.decode()))
 study=a.study
 for p in sorted((study/'inputs').rglob('*.json')):add('input-metadata',p,'study/'+str(p.relative_to(study)))
 for parent in ['native-complete','native']:
  for p in sorted((study/parent).rglob('*')):
   if p.is_file() and (p.suffix=='.log' or p.name in ['plan.json','receipt.json','commands.json','execution.json','prediction-freeze.json']):add('native-metadata',p,'study/'+str(p.relative_to(study)))
 for name in ['score','score-repeat','score-final','score-final-repeat','regression']:
  for p in sorted((study/name).rglob('*')):
   if p.is_file() and p.suffix in ['.json','.log']:add('checks-and-scores',p,'study/'+str(p.relative_to(study)))
 for p in sorted(study.iterdir()):
  if p.is_file() and (p.suffix in ['.log','.txt'] or p.name in ['score-repeat.json','score-final-repeat.json','external-artifacts.json','source-context.json']):add('execution-and-sources',p,'study/'+p.name)
 for p in sorted((study/'runtime').glob('*.sha256')):add('execution-and-sources',p,'study/runtime/'+p.name)
 for p in sorted(Path(__file__).parent.iterdir()):
  if p.is_file():add('executed-code',p,'repo/'+str(p.relative_to(a.repo)))
 for name in ['Sources/NumiVivoKit/Omics/VivoPerturbation.swift','Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift','Sources/NumiVivoCLI/VivoWorkflowCLIImplementation.swift','Tests/NumiVivoIntegrationTests/SingleCellPerturbationTests.swift','Tools/Omics/H5AD/build.sh','Tools/Omics/H5AD/CLIMain.swift','Tools/Omics/Reduction/GaussianKernel/TestMain.swift']:
  add('executed-code',a.repo/name,'repo/'+name)
 records=[]
 for group,files in sorted(groups.items()):
  b=(json.dumps(dict(schemaVersion=1,sourceFiles=files),sort_keys=True,separators=(',',':'))+'\n').encode();stored=gzip.compress(b,compresslevel=9,mtime=0);name=group+'.json.gz';(a.out/name).write_bytes(stored)
  records.append(dict(sourcePath=group+'.json',sourceBytes=len(b),sourceSHA256=digest(b),storedPath=name,storedBytes=len(stored),storedSHA256=digest(stored),gzipEncoded=True,bundledSourceFiles=True))
 (a.out/'manifest.json').write_text(json.dumps(dict(schemaVersion=1,records=records),indent=2)+'\n')
 print(json.dumps(dict(groups=len(records),bundledOriginalFiles=sum(len(v) for v in groups.values()),storedBytes=sum(x['storedBytes'] for x in records))))
if __name__=='__main__':main()
