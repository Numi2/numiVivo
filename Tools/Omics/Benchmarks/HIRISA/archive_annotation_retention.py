#!/usr/bin/env python3
"""Preserve all annotation diagnostic results, original attempts and source bindings."""
import argparse,gzip,hashlib
from pathlib import Path
from check_integration_response import sha
from prepare_annotation_retention import read,write


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('study','root','out'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();r=a.root;complete=read(r/'complete.json');freeze=read(r/'freeze.json')
 assert complete['allFivePerCellSVDChecksPassed'] and complete['allOriginalCellsScored'] and complete['identityExact']
 for name,h in complete['files'].items():assert sha(r/name)==h
 for name,h in freeze['files'].items():assert sha(r/name)==h
 for matrix in freeze['matrices'].values():assert sha(a.study/matrix['path'])==matrix['SHA256']
 unchanged=[]
 for name in ('baseline','identity','erasure','native','harmony-7','harmony-19','harmony-41'):
  old,new=read(r/'attempt-1'/(name+'.json')),read(r/(name+'.json'))
  assert len(old['folds'])==len(new['folds'])==114
  for x,y in zip(old['folds'],new['folds']):
   for k in ('id','confusion','metrics','accuracy','balancedRecall','trainingCounts','queryCounts'):assert x[k]==y[k],(name,x['id'],k)
  assert old['comparisons']==new['comparisons']
  unchanged.append(name)
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 paths=[('study/'+str(p.relative_to(r)),p) for p in r.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix in {'.py','.json','.log','.npz','.md'}]
 paths.append(('tools/archive_annotation_retention.py',Path(__file__)))
 for name,source in sorted(paths):
  assert not source.is_symlink() and source.stat().st_size<32*2**20
  raw=source.read_bytes();digest=hashlib.sha256(raw).hexdigest();encoded=gzip.compress(raw,compresslevel=6,mtime=0)
  target=a.out/(name+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(encoded);assert sha(source)==digest
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=digest,storedPath=name+'.gz',storedBytes=len(encoded),storedSHA256=sha(target),gzipEncoded=True))
 write(a.out/'manifest.json',dict(schemaVersion=1,records=records,cells=1612594,folds=114,labels=31,strata=23,completeSHA256=sha(r/'complete.json'),freezeSHA256=sha(r/'freeze.json'),externalMatrices=freeze['matrices'],externalSourceSHA256=freeze['sourceSHA256'],externalMetadataSHA256=freeze['metadataSHA256'],previousAttemptConfusionAndMetricsExact=unchanged,allFivePerCellSVDChecksPassed=True,scope='Complete original-cohort annotation retention diagnostic, all classes/folds and failed/insufficient gates. Author labels are nonauthoritative; no independent biological, prospective or general integration qualification. Original donor-subtraction attempt retained; explicit training-only summation leaves every confusion and metric unchanged.'))
 print(dict(members=len(records),storedBytes=sum(x['storedBytes'] for x in records),manifestSHA256=sha(a.out/'manifest.json')))
if __name__=='__main__':main()
