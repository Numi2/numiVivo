"""Advance verified folds through frozen prediction, reference and scoring recipes."""
from pathlib import Path
import fcntl,json,os,subprocess,sys,time
study,root=map(Path,sys.argv[1:3]);lock=(root/'advance.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);folds=json.loads((study/'protocol.json').read_text())['folds'];state={'status':'running','pid':os.getpid(),'startedUnix':time.time(),'completed':[]};env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
def read(p):return json.loads(p.read_text())
def save():
 t=root/'advance-state.tmp';t.write_text(json.dumps(state,indent=2)+'\n');t.replace(root/'advance-state.json')
try:
 while len(state['completed'])<len(folds):
  for fold in folds:
   tag=fold['tag'];out=root/tag
   if tag in state['completed']:continue
   watch=study/tag/fold['origin']/'verification-watch.json'
   if not watch.exists() or read(watch)['status']!='completed-all-workers-and-verification':continue
   p=out/'state.json'
   if p.exists() and read(p)['status']=='running':
    pid=read(p)['pid'];ps=subprocess.run(['ps','-p',str(pid),'-o','command='],capture_output=True,text=True)
    if ps.returncode==0 and str(root/'run.py') in ps.stdout:continue
    raise RuntimeError((tag,'prediction state running without its owner; inspect partial files'))
   state['activeFold']=tag;save()
   for recipe,result,expected in [('run.py','state.json','completed-control-only-predictions'),('verify.py','verification.json','passed-control-only-predictions'),('score.py','scores.json','scored-development-fold'),('verify_scores.py','scores-verification.json','passed-independent-sparse-scoring')]:
    target=out/result
    if target.exists() and read(target)['status']==expected:continue
    with (root/(tag+'-'+recipe+'.log')).open('ab') as log:subprocess.run([sys.executable,str(root/recipe),str(study),str(root),tag],stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
    assert read(target)['status']==expected
   state['completed'].append(tag);state.pop('activeFold',None);save();print('Scored and verified',tag,flush=True)
  if len(state['completed'])==len(folds):break
  parent=read(study/'fold-state.json');ps=subprocess.run(['ps','-p',str(parent['pid']),'-o','command='],capture_output=True,text=True)
  if parent['status']=='attention-required':raise RuntimeError('joint fitting scheduler needs attention; completed prediction folds retained')
  if parent['status']=='running' and not (ps.returncode==0 and 'run_folds.py' in ps.stdout and str(study) in ps.stdout):raise RuntimeError('joint scheduler owner absent; inspect fitting workers')
  if parent['status']=='completed-all-declared-folds':
   # An independently launched prediction may still own a fold.
   if not any((root/f['tag']/'state.json').exists() and read(root/f['tag']/'state.json')['status']=='running' for f in folds if f['tag'] not in state['completed']):raise RuntimeError('all fits terminal but an unadvanced fold remains')
  state['lastOwnerCheckUnix']=time.time();save();time.sleep(15)
 state['status']='completed-all-fold-prediction-and-scoring'
except BaseException as e:state.update(status='attention-required',error=repr(e));raise
finally:state['updatedUnix']=time.time();save()
