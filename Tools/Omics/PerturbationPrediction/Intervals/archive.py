#!/usr/bin/env python3
"""Retain exact compact interval evidence; full bundles are separately hash-bound."""
import argparse,gzip,hashlib,json
from pathlib import Path

def digest(b):return hashlib.sha256(b).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);groups={}
 def add(group,path,label):
  b=path.read_bytes();groups.setdefault(group,[]).append(dict(path=label,bytes=len(b),SHA256=digest(b),rawUTF8=b.decode()))
 for p in sorted((a.study/'native').rglob('*.json')):add('native-metadata',p,'study/'+str(p.relative_to(a.study)))
 for name in ['score','score-repeat','regression']:
  for p in sorted((a.study/name).rglob('*')):
   if p.is_file() and p.suffix in ['.json','.log']:add('checks-and-scores',p,'study/'+str(p.relative_to(a.study)))
 for p in sorted(a.study.iterdir()):
  if p.is_file() and p.suffix in ['.log','.json','.txt']:add('execution',p,'study/'+p.name)
 for p in sorted((a.study/'runtime').glob('*.sha256')):add('execution',p,'study/runtime/'+p.name)
 for p in sorted(Path(__file__).parent.iterdir()):
  if p.is_file():add('source-code',p,'repo/'+str(p.relative_to(a.repo)))
 for name in ['Sources/NumiVivoKit/Omics/VivoPerturbation.swift','Sources/NumiVivoKit/Omics/VivoOmicsLinearStatistics.swift','Tests/NumiVivoIntegrationTests/SingleCellPerturbationTests.swift','Tools/Omics/H5AD/build.sh','Tools/Omics/PerturbationPrediction/CrossStudyIFNB/test.sh','Tools/Omics/PerturbationPrediction/CrossStudyIFNB/regression.py','Tools/Omics/Reduction/GaussianKernel/TestMain.swift']:
  add('source-code',a.repo/name,'repo/'+name)
 records=[]
 for name,files in sorted(groups.items()):
  raw=(json.dumps(dict(schemaVersion=1,sourceFiles=files),sort_keys=True,separators=(',',':'))+'\n').encode();encoded=gzip.compress(raw,compresslevel=9,mtime=0);stored=name+'.json.gz';(a.out/stored).write_bytes(encoded)
  records.append(dict(sourcePath=name+'.json',sourceBytes=len(raw),sourceSHA256=digest(raw),storedPath=stored,storedBytes=len(encoded),storedSHA256=digest(encoded),gzipEncoded=True,bundledSourceFiles=True))
 (a.out/'manifest.json').write_text(json.dumps(dict(schemaVersion=1,records=records),indent=2)+'\n')
 print(json.dumps(dict(groups=len(groups),originalFiles=sum(len(v) for v in groups.values()),storedBytes=sum(x['storedBytes'] for x in records))))
if __name__=='__main__':main()
