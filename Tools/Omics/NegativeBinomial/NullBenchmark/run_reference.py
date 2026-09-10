#!/usr/bin/env python3
"""Run every frozen sham reference, retaining failed fits and compressed tables."""
import argparse,gzip,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--r-library',type=Path,required=True)
a=p.parse_args();owner=Path(__file__).resolve().parents[2]/'Bioconductor/run.R';records=[]
for study in ['kang','hagai']:
 for seed in range(1,11):
  d=a.root/study/str(seed);out=d/'r-reference';assert not out.exists()
  start=time.monotonic();cmd=['Rscript',str(owner),str(d/'r-input'),str(out)]
  with (d/'r-run.log').open('w') as log:
   r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'R_LIBS_USER':str(a.r_library),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
  records.append(dict(study=study,seed=seed,exitCode=r.returncode,seconds=time.monotonic()-start,command=cmd))
  for file in out.glob('*.tsv'):
   raw=file.read_bytes();compressed=gzip.compress(raw,mtime=0);dest=file.with_suffix('.tsv.gz');assert not dest.exists()
   dest.write_bytes(compressed);assert gzip.decompress(dest.read_bytes())==raw;file.unlink()
  (a.root/'reference-runs.json').write_text(json.dumps(dict(ownerSHA256=hashlib.sha256(owner.read_bytes()).hexdigest(),runs=records),indent=2)+'\n')
  print(json.dumps(records[-1]),flush=True)
