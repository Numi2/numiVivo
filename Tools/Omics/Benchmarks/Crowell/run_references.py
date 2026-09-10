#!/usr/bin/env python3
"""Run the pinned common-design Bioconductor comparison on every available population."""
import argparse,json,os,subprocess,time,hashlib
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','bioconductor','rscript','r_library','python']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();protocol=json.loads((a.root/'protocol.json').read_text());records=[]
env={**os.environ,'R_LIBS_USER':str(a.r_library),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
for case in protocol['cases']:
 d=a.root/case['id']
 if not (d/'reference-input').exists():continue
 for stage,command in [('reference',[a.rscript,a.bioconductor/'run.R',d/'reference-input',d/'reference']),('compare',[a.python,a.bioconductor/'compare.py','--input',d/'reference-input','--reference',d/'reference','--out',d/'comparison.json'])]:
  start=time.monotonic()
  with (d/(stage+'.log')).open('w') as f:r=subprocess.run(list(map(str,command)),stdout=f,stderr=subprocess.STDOUT,env=env)
  records.append(dict(case=case['id'],cellGroup=case['cellGroup'],stage=stage,exitCode=r.returncode,seconds=time.monotonic()-start,command=list(map(str,command))))
  (a.root/'reference-runs.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records[-1]),flush=True)
  if r.returncode:break
(a.root/'reference-complete.json').write_text(json.dumps(dict(status='all-available-populations-attempted',runs=records,failures=sum(r['exitCode']!=0 for r in records),runnerSHA256=hashlib.sha256((a.bioconductor/'run.R').read_bytes()).hexdigest(),comparisonSHA256=hashlib.sha256((a.bioconductor/'compare.py').read_bytes()).hexdigest()),indent=2)+'\n')
if any(r['exitCode'] for r in records):raise SystemExit(1)
