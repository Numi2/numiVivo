#!/usr/bin/env python3
"""Run native assay partition, full RNA aggregation and exact reconstruction."""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--binary',type=Path,required=True);a=p.parse_args();r=a.root;prepared=r/'prepared-complete';runs=[]
env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
assert json.loads((prepared/'inventory.json').read_text())['entries']==129839577
(r/'native-freeze.json').write_text(json.dumps(dict(binary=str(a.binary),binarySHA256=sha(a.binary),sourceSHA256=sha(prepared/'original.h5ad')),indent=2)+'\n')
def run(label,args):
 assert not (r/(label+'.log')).exists()
 with (r/(label+'.stdout')).open('wb') as out,(r/(label+'.log')).open('wb') as err:
  start=time.monotonic();p=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],env=env,stdout=out,stderr=err)
 runs.append(dict(label=label,arguments=list(map(str,args)),returnCode=p.returncode,seconds=time.monotonic()-start));(r/'native-progress.json').write_text(json.dumps(runs,indent=2)+'\n');assert p.returncode==0,(label,(r/(label+'.log')).read_text())
for assay in ['rna','guide']:
 run(assay+'-projection',['singlecell-h5ad-project',prepared/'original.h5ad','--plan',prepared/(assay+'-projection-plan.json'),'--output',r/(assay+'-projection')])
 run(assay+'-projection-verify',['singlecell-h5ad-project-verify',r/(assay+'-projection')])
run('aggregate',['singlecell-h5ad-pseudobulk',r/'rna-projection/projected.h5ad','--plan',prepared/'plan.json','--output',r/'aggregate'])
run('aggregate-verify',['singlecell-h5ad-pseudobulk-verify',r/'aggregate'])
(r/'native-checks.json').write_text(json.dumps(dict(status='passed',runs=runs,predictionFitted=False),indent=2)+'\n');print('All native assay partitions, full RNA aggregation and replays passed')
