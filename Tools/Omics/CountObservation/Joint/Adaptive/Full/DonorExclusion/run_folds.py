"""Run all declared folds with four native workers and independent checkers.

Frozen script copies prevent later repository edits from changing a live run.
Per-fold checkpoints and partial outputs are retained on any failure.
"""
from pathlib import Path
import collections,fcntl,hashlib,json,os,subprocess,sys,time
root,repo=map(Path,sys.argv[1:3]);lock=(root/'fold-controller.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
protocol=json.loads((root/'protocol.json').read_text());frozen=root/'runtime-recipes/Full';frozen.mkdir(parents=True,exist_ok=True);hashes={}
for name in ['run.py','verify.py','watch_verify.py']:
 p=frozen/name;source=repo/'Tools/Omics/CountObservation/Joint/Adaptive/Full'/name
 if not p.exists():p.write_bytes(source.read_bytes())
 hashes[name]=hashlib.sha256(p.read_bytes()).hexdigest()
freeze=root/'joint-runtime-recipes.json'
if freeze.exists():assert json.loads(freeze.read_text())==hashes
else:freeze.write_text(json.dumps(hashes,indent=2)+'\n')
folds=sorted(protocol['folds'],key=lambda f:(int(f['tag'].split('-')[1]),0 if f['origin']=='HIRISA' else 1));active={};done={};state={'status':'running','pid':os.getpid(),'startedUnix':time.time(),'maximumConcurrentNativeFits':4,'completed':{},'active':{}}
def save():
 state['active']={tag:{'fitPID':v['fit'].pid,'watchPID':v['watch'].pid if v.get('watch') else None,'root':v['fold']['root']} for tag,v in active.items()};state['completed']=done;state['updatedUnix']=time.time();p=root/'fold-state.tmp';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(root/'fold-state.json')
try:
 while folds or active:
  while folds and len(active)<4:
   disk=os.statvfs(root)
   if disk.f_bavail*disk.f_frsize<1024*1024*1024:raise OSError('insufficient free storage before starting another fold')
   fold=folds.pop(0);tag=fold['tag'];directory=Path(fold['root']);origin=fold['origin'];out=directory/origin;out.mkdir(exist_ok=True);existing=out/'verification-watch.json'
   if existing.exists() and json.loads(existing.read_text())['status']=='completed-all-workers-and-verification':
    v=json.loads((out/'verification-state.json').read_text());assert v['genesVerified']==fold['featureCount'] and not v['errors'];done[tag]={'status':'verified-complete','genes':v['genesVerified']};continue
   env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',NUMIVIVO_FULL_OUTPUT=str(out),NUMIVIVO_FULL_RUNTIME=str(directory/'runtime-freeze.json'),NUMIVIVO_FULL_SHARD_COUNT='1',NUMIVIVO_FULL_SHARD_INDEX='0')
   log=(directory/'joint-run.log').open('ab');fit=subprocess.Popen([sys.executable,str(frozen/'run.py'),str(directory),origin],stdout=log,stderr=subprocess.STDOUT,env=env);active[tag]={'fold':fold,'fit':fit,'fitLog':log,'env':env,'watch':None,'launched':time.monotonic()};save();print('Started',tag,'fit PID',fit.pid,flush=True)
  for tag,v in list(active.items()):
   fold=v['fold'];directory=Path(fold['root']);out=directory/fold['origin'];fit=v['fit'];state_file=out/'state.json'
   if v['watch'] is None and state_file.exists():
    st=json.loads(state_file.read_text())
    if st['pid']==fit.pid:
     log=(directory/'joint-watch.log').open('ab');v['watchLog']=log;v['watch']=subprocess.Popen([sys.executable,str(frozen/'watch_verify.py'),str(directory),fold['origin']],stdout=log,stderr=subprocess.STDOUT,env=v['env']);save()
   code=fit.poll()
   if code is not None and code!=0:raise RuntimeError((tag,'native controller failed',code))
   if v['watch'] is None:
    if time.monotonic()-v['launched']>60:raise RuntimeError((tag,'controller did not establish its state'))
    continue
   checked=v['watch'].poll()
   if checked is not None and checked!=0:raise RuntimeError((tag,'independent watcher failed',checked))
   if code==0 and checked==0:
    report=json.loads((out/'verification-state.json').read_text());watch=json.loads((out/'verification-watch.json').read_text());assert watch['status']=='completed-all-workers-and-verification' and report['genesVerified']==fold['featureCount'] and not report['errors'];done[tag]={'status':'verified-complete','genes':report['genesVerified'],'fitStates':report['fitStates']};v['fitLog'].close();v['watchLog'].close();del active[tag];save();print('Completed',tag,flush=True)
  save()
  if active:time.sleep(5)
 state['status']='completed-all-declared-folds'
except BaseException as e:state.update(status='attention-required',error=repr(e));raise
finally:state.update(finishedUnix=time.time(),seconds=time.time()-state['startedUnix']);save()
