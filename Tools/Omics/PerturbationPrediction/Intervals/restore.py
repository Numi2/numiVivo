#!/usr/bin/env python3
"""Restore one complete native interval bundle only after checking its frozen archive."""
import argparse,hashlib,json,os,tarfile,tempfile
from pathlib import Path,PurePosixPath

def main():
 p=argparse.ArgumentParser();p.add_argument('--native',type=Path,required=True);p.add_argument('--fold',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 freeze=json.loads((a.native/'prediction-freeze.json').read_text());matches=[f for f in freeze['folds'] if f['fold']==a.fold];assert len(matches)==1;record=matches[0]
 rel=PurePosixPath(record['path']);assert not rel.is_absolute() and '..' not in rel.parts
 archive=a.native/record['path'];assert archive.stat().st_size==record['bytes'] and hashlib.sha256(archive.read_bytes()).hexdigest()==record['SHA256']
 assert not a.out.exists();assert sum(x['bytes'] for x in record['members'].values())<=1_073_741_824
 staging=Path(tempfile.mkdtemp(prefix='.interval-restore-',dir=a.out.parent))
 with tarfile.open(archive,'r:gz') as t:
  members=t.getmembers();assert len(members)==len(record['members']) and {x.name for x in members}==set(record['members'])
  for member in members:
   name=PurePosixPath(member.name);expected=record['members'][member.name]
   assert member.isfile() and not name.is_absolute() and '..' not in name.parts and '\\' not in member.name and member.size==expected['bytes']
   target=staging.joinpath(*name.parts);target.parent.mkdir(parents=True,exist_ok=True);digest=hashlib.sha256();size=0
   with t.extractfile(member) as source,target.open('xb') as output:
    for chunk in iter(lambda:source.read(1_048_576),b''):output.write(chunk);digest.update(chunk);size+=len(chunk)
   assert size==expected['bytes'] and digest.hexdigest()==expected['SHA256'],member.name
 assert not a.out.exists();staging.rename(a.out)
 print(json.dumps(dict(status='restored',fold=a.fold,files=len(record['members']),output=str(a.out),archiveSHA256=record['SHA256'])))
if __name__=='__main__':main()
