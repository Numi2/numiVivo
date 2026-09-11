"""Predict a verified complete fold using only omitted-donor control counts."""
from pathlib import Path
import collections,fcntl,gzip,hashlib,json,os,shlex,subprocess,sys,time
import h5py,numpy as np
study,root=map(Path,sys.argv[1:3]);tag=sys.argv[3];root.mkdir(exist_ok=True);out=root/tag;out.mkdir(exist_ok=True)
lock=(out/'controller.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()
def encode(v):return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
fold=next(f for f in read(study/'protocol.json')['folds'] if f['tag']==tag);origin=fold['origin'];src=study/tag;fits=src/origin;m=read(src/(origin+'-manifest.json.gz'));report=read(fits/'verification-state.json');watch=read(fits/'verification-watch.json')
assert watch['status']=='completed-all-workers-and-verification' and report['genesVerified']==fold['featureCount'] and not report['errors']
assert m['excludedDonorID']==fold['excludedDonorID'] and sorted({g['donorID'] for g in m['groups']})==fold['trainingDonorIDs']
freeze=read(root/'runtime-freeze.json');cache=src/(origin+'-cells.h5');assert sha(cache)==m['cacheSHA256']
plan={'tag':tag,'excludedDonorID':m['excludedDonorID'],'modelManifestSHA256':sha(src/(origin+'-manifest.json.gz')),'cacheSHA256':m['cacheSHA256'],'binarySHA256':freeze['binarySHA256'],'plannedTreatedLibraryCounts':[10000],'treatedOutcomeInputs':False,'qualification':'Control-only plug-in joint-rate predictions. Standardized count scenario is not the observed treated library. No parameter uncertainty or calibrated intervals. Development donors only.'}
p=out/'protocol.json'
if p.exists():assert read(p)==plan
else:p.write_text(json.dumps(plan,indent=2)+'\n')
state={'status':'running','pid':os.getpid(),'startedUnix':time.time(),'genes':0,'states':{}}
def save():
 p=out/'state.tmp';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(out/'state.json')
counts=collections.Counter();save()
try:
 with h5py.File(cache,'r') as h:
  candidates=[k for k in h if h[k].attrs.get('donorID')==m['excludedDonorID'] and h[k].attrs.get('conditionID')=='control'];assert len(candidates)==1
  group=h[candidates[0]];libraries=group['depthLibraries'][:].tolist();cells=group['depthCells'][:].tolist();depth=group['depthIndex'][:];ptr=group['indptr'][:]
  for first in range(0,fold['featureCount'],64):
   last=min(first+64,fold['featureCount']);shard=f'{first:05d}-{last:05d}';r=read(fits/(shard+'-receipt.json'));model_path=fits/(shard+'-output.jsonl.gz');assert sha(model_path)==r['compressedOutputSHA256'];receipt=out/(shard+'-receipt.json');ip=out/(shard+'-input.jsonl.gz');op=out/(shard+'-output.jsonl.gz')
   if receipt.exists():
    rr=read(receipt);assert rr['modelSHA256']==sha(model_path) and rr['binarySHA256']==freeze['binarySHA256'] and rr['inputSHA256']==sha(ip) and rr['outputSHA256']==sha(op)
   else:
    assert os.statvfs(out).f_bavail*os.statvfs(out).f_frsize>512*1024*1024
    header={'origin':origin,'queryDonorID':m['excludedDonorID'],'trainingDonorIDs':fold['trainingDonorIDs'],'firstFeatureIndex':first,'featureIDs':[g['featureID'] for g in m['genes'][first:last]],'libraryCounts':libraries,'cellsPerLibrary':cells};lines=[encode(header)]
    with gzip.open(model_path,'rt') as f:
     for j,line in enumerate(f,first):
      g=json.loads(line);assert g['featureIndex']==j and g['featureID']==m['genes'][j]['featureID'];row={k:g[k] for k in ['featureIndex','featureID','status']};row.update(bins=[],counts=[])
      if g['status']=='boundedContinuousLikelihood':
       row['model']=g['model'];assert sorted(g['model']['model']['trainingDonorIDs'])==fold['trainingDonorIDs'];lo,hi=map(int,ptr[j:j+2]);v=np.bincount(depth[group['indices'][lo:hi]],weights=group['data'][lo:hi],minlength=len(libraries)).astype(np.uint64);bins=np.flatnonzero(v);row.update(bins=bins.tolist(),counts=v[bins].tolist())
      lines.append(encode(row))
    assert len(lines)==last-first+1;data=b'\n'.join(lines)+b'\n';packed=gzip.compress(data,mtime=0)
    if ip.exists():assert ip.read_bytes()==packed
    else:ip.write_bytes(packed)
    native=freeze['nativePath'];check='import hashlib,pathlib;assert hashlib.sha256(pathlib.Path('+repr(native)+').read_bytes()).hexdigest()=='+repr(freeze['binarySHA256']);command='python3 -c '+shlex.quote(check)+' && '+shlex.quote(native);t=time.time()
    result=subprocess.run(['ssh','-4','-C','-o','ConnectTimeout=10','-o','ServerAliveInterval=10','-o','ServerAliveCountMax=1','macmini',command],input=data,capture_output=True)
    (out/(shard+'-stderr.log')).write_bytes(result.stderr)
    if result.returncode:(out/(shard+'-failed-output')).write_bytes(result.stdout);result.check_returncode()
    rows=[json.loads(v) for v in result.stdout.splitlines()];assert len(rows)==last-first
    for j,g in enumerate(rows,first):assert g['featureIndex']==j and g['featureID']==m['genes'][j]['featureID']
    assert not op.exists();op.write_bytes(gzip.compress(result.stdout,mtime=0));receipt.write_text(json.dumps({'status':'completed','modelSHA256':sha(model_path),'binarySHA256':freeze['binarySHA256'],'inputSHA256':sha(ip),'outputSHA256':sha(op),'seconds':time.time()-t},indent=2)+'\n')
   with gzip.open(op,'rt') as f:
    for line in f:counts[json.loads(line)['status']]+=1
   state.update(genes=last,states=dict(counts),lastShard=shard);save()
 state['status']='completed-control-only-predictions'
except BaseException as e:state.update(status='failed',error=repr(e));raise
finally:state.update(finishedUnix=time.time());save()
