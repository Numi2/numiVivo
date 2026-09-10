import os,json,time,subprocess
from pathlib import Path
r=Path('/Users/home/numivivo-empirical-prior-20260910');remote='/Users/n/numivivo-empirical-prior-20260910';repo='/Users/home/numivivo-single-cell-20260909';python='/Users/home/numivivo-singlecell-benchmark-py/bin/python'
records=json.loads((r/'protocol.json').read_text())['records'];env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
for rec in records:
 name=rec['id'];assert all(c.isalnum() or c=='-' for c in name);start=time.monotonic()
 while True:
  t=subprocess.run(['ssh','macmini',f'test -f {remote}/{name}/native/archive.json && test ! -f {remote}/{name}/native/report.json'])
  if t.returncode==0:break
  if time.monotonic()-start>1200:raise RuntimeError('Archive observation timed out; native job state must be inspected before further action')
  time.sleep(5)
 subprocess.run(['rsync','-a',f'macmini:{remote}/{name}/native',str(r/name)+'/'],check=True)
 with (r/(name+'-check.log')).open('w') as log:
  done=subprocess.run([python,'Tools/Omics/NegativeBinomial/EmpiricalPrior/check.py','--root',str(r),'--case',name],cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT)
 print(json.dumps(dict(case=name,checkerExit=done.returncode)),flush=True)
 if done.returncode:raise RuntimeError('Retain failed check and inspect before proceeding')
results=[json.loads((r/('check-'+x['id']+'.json')).read_text()) for x in records]
assert len({x['checkerSHA256'] for x in results})==1
combined={**results[0],'studies':[s for x in results for s in x['studies']],'seconds':sum(x['seconds'] for x in results),'failedFeatures':sum(x['failedFeatures'] for x in results)}
(r/'check.json').write_text(json.dumps(combined,sort_keys=True,indent=2)+'\n')
