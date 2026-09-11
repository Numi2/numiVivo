"""Reconstruct query inputs from original control cells and verify native outputs."""
from pathlib import Path
import collections,gzip,hashlib,json,sys,time
import h5py,numpy as np
from reference import likelihood,moments,sampling
study,root=map(Path,sys.argv[1:3]);tag=sys.argv[3];out=root/tag
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()
fold=next(f for f in read(study/'protocol.json')['folds'] if f['tag']==tag);origin=fold['origin'];src=study/tag;manifest=read(src/(origin+'-manifest.json.gz'));plan=read(out/'protocol.json');state=read(out/'state.json');assert state['status']=='completed-control-only-predictions';assert plan['treatedOutcomeInputs'] is False and plan['plannedTreatedLibraryCounts']==[10000]
cache=src/(origin+'-cells.h5');assert sha(cache)==manifest['cacheSHA256']==plan['cacheSHA256'];assert sha(src/(origin+'-manifest.json.gz'))==plan['modelManifestSHA256'];errors=[];maximum=0.;compared=0;statuses=collections.Counter();degenerate=underflowed=0;bindings={};start=time.time()
def close(a,b,label):
 global maximum,compared
 a=np.asarray(a,float);b=np.asarray(b,float);assert a.shape==b.shape;err=float(np.max(np.abs(a-b)/np.maximum(1,np.abs(b)))) if a.size else 0.;maximum=max(maximum,err);compared+=int(a.size)
 if not np.isfinite(err) or err>2e-7:errors.append({'where':label,'maximumScaledError':err})
with h5py.File(cache,'r') as h:
 keys=[k for k in h if h[k].attrs.get('donorID')==fold['excludedDonorID'] and h[k].attrs.get('conditionID')=='control'];assert len(keys)==1;group=h[keys[0]];libraries=group['libraryCounts'][:];unique,inverse,n=np.unique(libraries,return_inverse=True,return_counts=True);ptr=group['indptr'][:]
 for first in range(0,fold['featureCount'],64):
  last=min(first+64,fold['featureCount']);shard=f'{first:05d}-{last:05d}';ip=out/(shard+'-input.jsonl.gz');op=out/(shard+'-output.jsonl.gz');receipt=read(out/(shard+'-receipt.json'));modelpath=src/origin/(shard+'-output.jsonl.gz');fitreceipt=read(src/origin/(shard+'-receipt.json'));assert sha(modelpath)==receipt['modelSHA256']==fitreceipt['compressedOutputSHA256'];assert sha(ip)==receipt['inputSHA256'] and sha(op)==receipt['outputSHA256'];assert receipt['binarySHA256']==plan['binarySHA256'];bindings[shard]={'inputSHA256':sha(ip),'outputSHA256':sha(op),'modelSHA256':sha(modelpath)}
  lines=gzip.decompress(ip.read_bytes()).splitlines();head=json.loads(lines[0]);inputs=[json.loads(l) for l in lines[1:]];outputs=[json.loads(l) for l in gzip.decompress(op.read_bytes()).splitlines()];models=[json.loads(l) for l in gzip.decompress(modelpath.read_bytes()).splitlines()];assert len(inputs)==len(outputs)==len(models)==last-first
  assert head['queryDonorID']==fold['excludedDonorID'] and head['trainingDonorIDs']==fold['trainingDonorIDs'];assert head['libraryCounts']==unique.tolist() and head['cellsPerLibrary']==n.tolist();assert head['firstFeatureIndex']==first and head['featureIDs']==[g['featureID'] for g in manifest['genes'][first:last]]
  for offset,(inp,row,fit) in enumerate(zip(inputs,outputs,models)):
   j=first+offset;feature=manifest['genes'][j]['featureID'];assert inp['featureID']==row['featureID']==fit['featureID']==feature;assert inp['featureIndex']==row['featureIndex']==fit['featureIndex']==j;assert inp['status']==fit['status'];statuses[row['status']]+=1
   if fit['status']!='boundedContinuousLikelihood':assert row['status']==fit['status'] and 'prediction' not in row and inp['bins']==inp['counts']==[];continue
   assert inp['model']==fit['model'];model=fit['model']['model'];assert model['trainingDonorIDs']==fold['trainingDonorIDs'] and fold['excludedDonorID'] not in model['trainingDonorIDs'];y=np.zeros(len(libraries),dtype=np.uint64);lo,hi=map(int,ptr[j:j+2]);y[group['indices'][lo:hi]]=group['data'][lo:hi];totals=np.bincount(inverse,weights=y).astype(np.uint64);bins=np.flatnonzero(totals);assert inp['bins']==bins.tolist() and inp['counts']==totals[bins].tolist()
   p=row['prediction'];assert row['status']=='conditionalPrediction' and p['featureID']==feature and p['queryDonorID']==fold['excludedDonorID'] and p['treatedConditionID']=='IFNB';assert bytes(p['querySource']['bytes']).hex()==hashlib.sha256(lines[0]+lines[offset+1]).hexdigest()
   _,ell,_=likelihood(y,libraries,model['controlCellDispersion']);cr=np.array(model['controlRatesCPM']);tr=np.array(model['treatedRatesCPM']);w=np.array(model['probabilities']);logs=np.full(w.shape,-np.inf);positive=w>0;logs[positive]=np.log(w[positive])+np.repeat([ell(r) for r in cr],len(tr))[positive];ww=np.exp(logs-np.max(logs));ww/=ww.sum();close(p['probabilities'],ww,feature+'/posterior');mm=moments(cr,tr,ww)
   for k,v in mm.items():close(p['moments'][k],v,feature+'/'+k)
   for k,v in sampling(mm,model['treatedCellDispersion'],[10000]).items():close(p['plannedTreatedCountMoments'][k],v,feature+'/'+k)
   support=np.isfinite(logs);rates=np.tile(tr,len(cr))[support];point=bool(np.all(rates==rates[0]));assert p['degenerateTreatedRateDistribution']==point;degenerate+=int(point);underflowed+=p['underflowedPosteriorComponents'];assert p['underflowedPosteriorComponents']==int(np.sum(support & (ww==0)))
report={'status':'passed-control-only-predictions' if not errors else 'failed','tag':tag,'genes':sum(statuses.values()),'states':dict(statuses),'numericalValuesCompared':compared,'maximumScaledError':maximum,'pointMassPredictions':degenerate,'underflowedPosteriorComponents':underflowed,'errors':errors,'shardBindings':bindings,'seconds':time.time()-start,'treatedOutcomeReads':False};tmp=out/'verification.tmp';tmp.write_text(json.dumps(report,indent=2)+'\n');tmp.replace(out/'verification.json');print({k:v for k,v in report.items() if k!='shardBindings'});assert not errors
