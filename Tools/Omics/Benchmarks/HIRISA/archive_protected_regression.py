#!/usr/bin/env python3
"""Archive the complete small evidence and bind externally retained coordinate files."""
import argparse,gzip,hashlib
from pathlib import Path
from check_integration_response import sha
from prepare_annotation_retention import read,write


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('study','root','out'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();r=a.root;c=read(r/'complete.json');o=read(r/'output-freeze.json');f=read(r/'execution-freeze.json');summary=read(r/'summary.json');owner=Path(__file__).parent
 assert c['allNumericalChecksPassed'] and c['allOriginalCellsPreserved']
 assert sha(r/'output-freeze.json')==c['outputFreezeSHA256'] and sha(r/'complete.json')==summary['completeSHA256']
 for n,h in c['results'].items():assert sha(r/'evaluation'/n)==h
 for n,h in o['files'].items():assert sha(r/n)==h
 for v in f['inputs'].values():assert sha(a.study/v['path'])==v['SHA256']
 execution=read(r/'evaluation-execution-freeze.json')
 # Preserve every Python source bound at evaluation, including imported decoders.
 for n,h in execution['sourceFiles'].items():assert sha(owner/n)==h
 a.out.mkdir(parents=True,exist_ok=False)
 paths=[('study/'+str(p.relative_to(r)),p) for p in r.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix in {'.json','.md','.log','.npz','.py'}]
 paths += [('tools/'+n,owner/n) for n in sorted(set(execution['sourceFiles'])|{'summarize_protected_regression.py','archive_protected_regression.py'})]
 records=[]
 for n,p in sorted(paths):
  assert not p.is_symlink() and p.stat().st_size<64*2**20
  raw=p.read_bytes();h=hashlib.sha256(raw).hexdigest();encoded=gzip.compress(raw,compresslevel=6,mtime=0)
  q=a.out/(n+'.gz');q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(encoded);assert sha(p)==h
  records.append(dict(sourcePath=n,sourceBytes=len(raw),sourceSHA256=h,storedPath=n+'.gz',storedBytes=len(encoded),storedSHA256=sha(q),gzipEncoded=True))
 write(a.out/'manifest.json',dict(schemaVersion=1,records=records,completeSHA256=sha(r/'complete.json'),summarySHA256=sha(r/'summary.json'),externalInputs=f['inputs'],externalOutputs={n:dict(path=str((r/n).relative_to(a.study)),bytes=(r/n).stat().st_size,SHA256=h) for n,h in o['files'].items() if n.endswith('.bin')},scope='All models, frozen protocol, original-row metadata, evaluator sources, results and independent checks archived. Original source and full coordinate inputs/outputs remain external and hash-bound; restore under the recorded study-relative paths. Failed biological preservation retained; no native promotion.'))
 print(dict(members=len(records),storedBytes=sum(v['storedBytes'] for v in records),manifestSHA256=sha(a.out/'manifest.json')))
if __name__=='__main__':main()
