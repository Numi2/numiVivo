#!/usr/bin/env python3
"""Run frozen full-cohort annotation diagnostics and every-fold SVD checks."""
import argparse,os,subprocess,time
from pathlib import Path
from prepare_annotation_retention import read,write
from check_integration_response import sha


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('study','root','python'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();r=a.root;owner=Path(__file__).parent
 freeze=read(r/'freeze.json');assert freeze['cells']==1612594
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
 commands=[]
 def run(name,args):
  start=time.time();command=[str(a.python),*map(str,args)]
  with (r/(name+'.log')).open('x') as stream:
   process=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,env=env)
   write(r/(name+'-process.json'),dict(pid=process.pid,driverPID=os.getpid(),command=command,startedUnix=start))
   code=process.wait()
  record=dict(name=name,returnCode=code,seconds=time.time()-start);commands.append(record);write(r/(name+'-status.json'),record);print(record,flush=True);assert code==0,name
 run('tests',[owner/'test_annotation_retention.py'])
 write(r/'execution-freeze.json',dict(startedUnix=time.time(),pid=os.getpid(),freezeSHA256=sha(r/'freeze.json'),sourceFiles={p.name:sha(p) for p in owner.glob('*annotation_retention*.py')},readerSHA256=sha(owner/'check_integration_response.py'),testsSHA256=sha(r/'tests.log'),scope='No integration refit. Fixed annotation retention evaluation of all original cells on five qualified matrices; author labels are nonauthoritative.'))
 for name,matrix,erase in [('baseline','baseline',False),('identity','baseline',False),('erasure','baseline',True)]+[(s,s,False) for s in ('native','harmony-7','harmony-19','harmony-41')]:
  args=[owner/'check_annotation_retention.py','--study',a.study,'--root',r,'--out',r/(name+'.json'),'--matrix',matrix]
  if name!='baseline':args+=['--baseline',r/'baseline.json']
  if erase:args+=['--erase']
  run(name,args)
 assert read(r/'baseline.json')==read(r/'identity.json'),'Exact identity reconstruction'
 for name in ('baseline','native','harmony-7','harmony-19','harmony-41'):
  run(name+'-oracle',[owner/'check_annotation_retention_oracle.py','--study',a.study,'--root',r,'--result',r/(name+'.json'),'--out',r/(name+'-oracle.json')])
 write(r/'complete.json',dict(status='completed-diagnostic',commands=commands,identityExact=True,allOriginalCellsScored=True,allFivePerCellSVDChecksPassed=True,freezeSHA256=sha(r/'freeze.json'),files={n:sha(r/n) for n in [*(s+'.json' for s in ('baseline','identity','erasure','native','harmony-7','harmony-19','harmony-41')),*(s+'-oracle.json' for s in ('baseline','native','harmony-7','harmony-19','harmony-41'))]}))
if __name__=='__main__':main()
