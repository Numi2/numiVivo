#!/usr/bin/env python3
"""Run the declared optimization sensitivity while retaining original references."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','null-root','r-library']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();owner=Path(__file__).resolve().parent;output=a.root/'sensitivity';output.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
(output/'declaration.json').write_text(json.dumps({name:sha(owner/name) for name in ['SENSITIVITY_PROTOCOL.md','sensitivity.R','sensitivity.py']},indent=2)+'\n')
records=[]
for study in ['kang','hagai']:
 for seed in range(1,11):
  d=a.null_root/study/str(seed);audit=a.root/study/str(seed);out=output/study/str(seed);out.parent.mkdir(parents=True,exist_ok=True)
  with (out.parent/(str(seed)+'.log')).open('w') as log:
   r=subprocess.run(['/opt/homebrew/bin/Rscript',str(owner/'sensitivity.R'),str(d/'r-input'),str(audit),str(d/'r-reference/native_size_factors-DESeq2.tsv.gz'),str(out)],stdout=log,stderr=subprocess.STDOUT,
    env={**os.environ,'R_LIBS_USER':str(a.r_library),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
  record=dict(study=study,seed=seed,exitCode=r.returncode,inputCheckSHA256=sha(audit/'checked.tsv.gz'))
  if (out/'status.json').exists():record['result']=json.loads((out/'status.json').read_text())
  records.append(record);(output/'runs.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
  print(json.dumps(dict(study=study,seed=seed,exitCode=r.returncode,result=record.get('result'))),flush=True)
assert len(records)==20 and all(r['exitCode']==0 and r['result']['status']=='completed' for r in records)
(output/'complete.json').write_text(json.dumps(dict(status='completed-all-twenty-optimization-sensitivities',baselineAndCandidateFits=40),indent=2)+'\n')
