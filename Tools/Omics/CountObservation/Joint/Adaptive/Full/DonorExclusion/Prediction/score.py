"""Score only complete, independently checked prediction folds; no fitting here."""
from pathlib import Path
import gzip,hashlib,json,sys,time
import h5py,numpy as np
study,root=map(Path,sys.argv[1:3]);tag=sys.argv[3];out=root/tag

def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()
fold=next(f for f in read(study/'protocol.json')['folds'] if f['tag']==tag);origin=fold['origin'];src=study/tag;manifest=read(src/(origin+'-manifest.json.gz'));v=read(out/'verification.json');assert v['status']=='passed-control-only-predictions' and v['genes']==fold['featureCount'] and not v['errors'];assert not (out/'scores.json').exists()
# Verify all prediction bytes before opening any outcome-bearing count array.
for shard,b in v['shardBindings'].items():
 assert sha(out/(shard+'-input.jsonl.gz'))==b['inputSHA256'] and sha(out/(shard+'-output.jsonl.gz'))==b['outputSHA256']
plan=read(root/'scoring-plan.json');cache=src/(origin+'-cells.h5');assert sha(cache)==manifest['cacheSHA256'];n=fold['featureCount'];pseudobulk={}
with h5py.File(cache,'r') as h:
 for key in h:
  if not isinstance(h[key],h5py.Group):continue
  g=h[key];donor=str(g.attrs['donorID']);condition=str(g.attrs['conditionID']);ptr=g['indptr'][:];totals=np.zeros(n,dtype=np.uint64);libraries=int(g['libraryCounts'][:].sum())
  for first in range(0,n,64):
   last=min(first+64,n);lo,hi=int(ptr[first]),int(ptr[last]);data=g['data'][lo:hi];s=np.concatenate([np.zeros(1,dtype=np.uint64),np.cumsum(data,dtype=np.uint64)]);totals[first:last]=s[ptr[first+1:last+1]-lo]-s[ptr[first:last]-lo]
  pseudobulk[(donor,condition)]=np.log1p(totals.astype(float)*1e6/libraries)
query=fold['excludedDonorID'];control=pseudobulk[(query,'control')];observed=pseudobulk[(query,'IFNB')];training=fold['trainingDonorIDs'];baseline=np.maximum(0,control+np.mean([pseudobulk[(d,'IFNB')]-pseudobulk[(d,'control')] for d in training],axis=0));prediction=np.full(n,np.nan);statuses=[];point=0
for shard in sorted(v['shardBindings']):
 with gzip.open(out/(shard+'-output.jsonl.gz'),'rt') as f:
  for line in f:
   row=json.loads(line);j=row['featureIndex'];assert j==len(statuses) and row['featureID']==manifest['genes'][j]['featureID'];statuses.append(row['status'])
   if row['status']=='conditionalPrediction':prediction[j]=np.log1p(row['prediction']['moments']['treatedMeanCPM']);point+=int(row['prediction']['degenerateTreatedRateDistribution'])
assert len(statuses)==n;mask=np.isfinite(prediction);assert int(mask.sum())==v['states']['conditionalPrediction']
def metrics(x):
 e=x[mask]-observed[mask];return {'genes':int(mask.sum()),'RMSE':float(np.sqrt(np.mean(e*e))),'MAE':float(np.mean(abs(e))),'sumSquaredError':float(e@e)}
models={name:metrics(x) for name,x in [('joint',prediction),('noChange',control),('trainingMeanResponse',baseline)]};gains={name:1-models['joint']['RMSE']/models[name]['RMSE'] for name in ['noChange','trainingMeanResponse']}
rows=[]
for j,status in enumerate(statuses):rows.append({'featureIndex':j,'featureID':manifest['genes'][j]['featureID'],'status':status,'observedLog1pCPM':float(observed[j]),'controlLog1pCPM':float(control[j]),'trainingMeanResponseLog1pCPM':float(baseline[j]),'predictionLog1pCPM':float(prediction[j]) if mask[j] else None})
raw=json.dumps(rows,separators=(',',':'),allow_nan=False).encode();table=out/'scores-by-gene.json.gz';table.write_bytes(gzip.compress(raw,mtime=0));report={'tag':tag,'status':'scored-development-fold','queryDonorID':query,'sourceGenes':n,'scoredGenes':int(mask.sum()),'excludedStates':{k:c for k,c in v['states'].items() if k!='conditionalPrediction'},'metrics':models,'relativeRMSEGain':gains,'fractionGenesBetterThanNoChange':float(np.mean(abs(prediction[mask]-observed[mask])<abs(control[mask]-observed[mask]))),'pointMassPredictions':point,'scoresSHA256':sha(table),'verifiedPredictionManifestSHA256':sha(out/'verification.json'),'scoringPlanSHA256':sha(root/'scoring-plan.json'),'finishedUnix':time.time(),'qualification':'One reused development donor; no all-fold success or calibrated interval claim.'};(out/'scores.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
