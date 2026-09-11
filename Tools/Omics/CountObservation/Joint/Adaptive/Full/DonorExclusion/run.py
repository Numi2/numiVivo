"""Execute every frozen training-only calibration on the pinned native owner."""
from pathlib import Path
import collections,gzip,hashlib,json,os,shlex,subprocess,sys,time
root=Path(sys.argv[1]);protocol=json.loads((root/'protocol.json').read_text());runtime=json.loads((root/'runtime-freeze.json').read_text());out=root/'native';out.mkdir(exist_ok=True)
state={'status':'running','pid':os.getpid(),'startedUnix':time.time(),'foldsCompleted':0,'activeFold':None}
def save():
 p=root/'state.tmp';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(root/'state.json')
save()
try:
 for fold in protocol['folds']:
  tag=fold['tag'];state['activeFold']=tag;save();p=out/(tag+'.json.gz');assert not p.exists();raw=gzip.decompress((root/'inputs'/(tag+'.json.gz')).read_bytes());assert hashlib.sha256(raw).hexdigest()==fold['inputSHA256'];partial=p.with_suffix('.partial.gz');assert not partial.exists()
  verify='import hashlib,pathlib;assert hashlib.sha256(pathlib.Path('+repr(runtime['nativePath'])+').read_bytes()).hexdigest()=='+repr(runtime['binarySHA256'])
  command='python3 -c '+shlex.quote(verify)+' && bash -o pipefail -c '+shlex.quote('/usr/bin/time -l '+shlex.quote(runtime['nativePath'])+' | gzip -c')
  t=time.time();r=subprocess.run(['ssh','-4','-C','-o','ConnectTimeout=10','-o','ServerAliveInterval=10','-o','ServerAliveCountMax=1','macmini',command],input=raw,capture_output=True);(out/(tag+'-stderr.log')).write_bytes(r.stderr);partial.write_bytes(r.stdout);r.check_returncode();data=gzip.decompress(r.stdout);model=json.loads(data);assert model['inputSHA256']==fold['inputSHA256'] and model['excludedDonorID']==fold['excludedDonorID']
  for side in ['control','treated']:assert model[side]['trainingDonorIDs']==fold['trainingDonorIDs'] and len(model[side]['features'])==fold['featureCount']
  partial.rename(p);receipt=dict(fold,status='completed',binarySHA256=runtime['binarySHA256'],outputSHA256=hashlib.sha256(data).hexdigest(),compressedOutputSHA256=hashlib.sha256(r.stdout).hexdigest(),seconds=time.time()-t,states={side:dict(collections.Counter(f['status'] for f in model[side]['features'])) for side in ['control','treated']});(out/(tag+'-receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n');state['foldsCompleted']+=1;save();print(json.dumps({'tag':tag,'seconds':receipt['seconds'],'states':receipt['states']}),flush=True)
 state['status']='completed-all-training-only-folds';state['activeFold']=None
except BaseException as e:state.update(status='failed',error=repr(e));raise
finally:state.update(finishedUnix=time.time(),seconds=time.time()-state['startedUnix']);save()
