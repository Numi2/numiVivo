#!/usr/bin/env python3
"""Freeze the four existing complete donor splits for native logistic qualification."""
import argparse,gzip,hashlib,json,shutil
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--prior',type=Path,required=True);p.add_argument('--native-prior',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 links=json.loads((a.repo/'Tools/Omics/ReferenceMapping/evidence/2026-09-09-native/source-links.json').read_text());sources={}
 for donor in ['human1','human2','human3','human4']:
  root=a.out/donor;root.mkdir();fit=json.loads((a.native_prior/donor/'fit.json').read_text());fit['logistic']=dict(penalty=1,gradientTolerance=1e-7,maximumIterations=2000,maximumWork=20000000000);write(root/'fit.json',fit)
  shutil.copy2(a.native_prior/donor/'query.json',root/'query.json')
  for group in ['train','query']:
   source=a.prior/donor/group/'projected.h5ad';h=sha(source);assert h==links[donor+'/'+group]['sourceSHA256']
   target=root/(group+'.h5ad.gz')
   with source.open('rb') as f,target.open('wb') as g:
    with gzip.GzipFile(fileobj=g,mode='wb',mtime=0,compresslevel=6) as z:shutil.copyfileobj(f,z,1048576)
   with gzip.open(target,'rb') as f:
    v=hashlib.sha256()
    for b in iter(lambda:f.read(1048576),b''):v.update(b)
   assert v.hexdigest()==h
   sources[donor+'/'+group]=dict(sourcePath=str(source),decodedSHA256=h,decodedBytes=source.stat().st_size,gzipPath=str(target.relative_to(a.out)),gzipSHA256=sha(target),gzipBytes=target.stat().st_size,priorSourceLink=links[donor+'/'+group])
  print(json.dumps(dict(prepared=donor)),flush=True)
 write(a.out/'input-freeze.json',dict(schemaVersion=1,protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),preparerSHA256=sha(__file__),sources=sources,files={str(p.relative_to(a.out)):sha(p) for p in sorted(a.out.rglob('*')) if p.is_file()},allFourDonors=True,fitStarted=False,queryLabelsExcludedByMapping=True))
if __name__=='__main__':main()
