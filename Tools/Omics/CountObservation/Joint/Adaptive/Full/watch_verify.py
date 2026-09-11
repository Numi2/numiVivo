"""Verify newly completed native shards while their owned workers continue."""
from pathlib import Path
import fcntl,json,os,subprocess,sys,time
root=Path(sys.argv[1]);origin=sys.argv[2];folder=Path(os.environ.get('NUMIVIVO_FULL_OUTPUT',str(root/origin)));worker_count=int(os.environ.get('NUMIVIVO_FULL_SHARD_COUNT','1'))
lock=(folder/'verification-watch.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
state={'status':'running','pid':os.getpid(),'origin':origin,'startedUnix':time.time(),'verifiedReceipts':0,'workerCount':worker_count}
def save():(folder/'verification-watch.json').write_text(json.dumps(state,indent=2)+'\n')
save();seen=-1
try:
 while True:
  count=len(list(folder.glob('*-receipt.json')))
  worker_states=[];live=[]
  for i in range(worker_count):
   p=folder/('state.json' if worker_count==1 else f'state-worker-{i}.json')
   if not p.exists():continue
   s=json.loads(p.read_bytes());worker_states.append(s)
   result=subprocess.run(['ps','-p',str(s['pid']),'-o','args='],capture_output=True,text=True)
   if result.returncode==0 and '/Full/run.py' in result.stdout and str(root) in result.stdout and origin in result.stdout:live.append(s['pid'])
  if count>seen and count>0:
   with (folder/'verification-watch.log').open('a') as log:
    p=subprocess.Popen([sys.executable,str(Path(__file__).with_name('verify.py')),str(root),origin],stdout=log,stderr=subprocess.STDOUT)
    state['activeVerifierPID']=p.pid;state['liveWorkerPIDs']=live;save();code=p.wait();state['activeVerifierPID']=None
   if code:state.update(status='failed-independent-verification',exitCode=code);break
   verified=json.loads((folder/'verification-state.json').read_bytes());assert verified['status']=='verified-completed-shards' and not verified['errors'];seen=count;state.update(verifiedReceipts=count,genesVerified=verified['genesVerified'],numericalValuesCompared=verified['numericalValuesCompared'],lastUpdateUnix=time.time());save();print(json.dumps(state),flush=True)
  completed=len(worker_states)==worker_count and all(s['status'] in ['completed-all-source-genes','completed-assigned-source-shards'] for s in worker_states)
  if not live:
   if completed and count==len(list(folder.glob('*-receipt.json'))) and count==seen:state['status']='completed-all-workers-and-verification';break
   state.update(status='worker-recovery-required',workerStates=[s['status'] for s in worker_states]);break
  time.sleep(10)
finally:state['finishedUnix']=time.time();state['liveWorkerPIDs']=live;save()
sys.exit(0 if state['status']=='completed-all-workers-and-verification' else 1)
