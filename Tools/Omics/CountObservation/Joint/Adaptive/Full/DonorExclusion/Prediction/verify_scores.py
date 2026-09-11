"""Check every scored value with SciPy CSC sums and scalar metric reductions."""
from pathlib import Path
import gzip,hashlib,json,math,sys
import h5py,numpy as np
from scipy import sparse
study,root=map(Path,sys.argv[1:3]);tag=sys.argv[3];out=root/tag
fold=next(f for f in json.loads((study/'protocol.json').read_text())['folds'] if f['tag']==tag);report=json.loads((out/'scores.json').read_text());table=out/'scores-by-gene.json.gz';assert hashlib.sha256(table.read_bytes()).hexdigest()==report['scoresSHA256'];rows=json.loads(gzip.decompress(table.read_bytes()));n=fold['featureCount'];assert len(rows)==n;raw={};compared=0;maximum=0.
def close(a,b):
 global compared,maximum
 err=abs(a-b)/max(1,abs(b));maximum=max(maximum,err);compared+=1;assert math.isfinite(err) and err<1e-12,(a,b,err)
with h5py.File(study/tag/(fold['origin']+'-cells.h5'),'r') as h:
 for key in h:
  g=h[key]
  if not isinstance(g,h5py.Group):continue
  ptr=g['indptr'][:];libs=g['libraryCounts'][:];total=sum(map(int,libs));values=[]
  for first in range(0,n,64):
   last=min(first+64,n);lo,hi=map(int,[ptr[first],ptr[last]]);x=sparse.csc_matrix((g['data'][lo:hi],g['indices'][lo:hi],ptr[first:last+1]-lo),shape=(len(libs),last-first));sums=np.asarray(x.sum(axis=0)).ravel();values.extend(math.log1p(int(v)*1e6/total) for v in sums)
  raw[(str(g.attrs['donorID']),str(g.attrs['conditionID']))]=values
query=fold['excludedDonorID'];donors=fold['trainingDonorIDs'];predicted={}
for p in sorted(out.glob('*-output.jsonl.gz')):
 for line in gzip.decompress(p.read_bytes()).splitlines():
  v=json.loads(line)
  if v['status']=='conditionalPrediction':predicted[v['featureIndex']]=math.log1p(v['prediction']['moments']['treatedMeanCPM'])
errors={k:[] for k in report['metrics']}
for j,row in enumerate(rows):
 y=raw[(query,'IFNB')][j];c=raw[(query,'control')][j];b=max(0,c+math.fsum(raw[(d,'IFNB')][j]-raw[(d,'control')][j] for d in donors)/len(donors));close(row['observedLog1pCPM'],y);close(row['controlLog1pCPM'],c);close(row['trainingMeanResponseLog1pCPM'],b)
 if j not in predicted:assert row['predictionLog1pCPM'] is None;continue
 close(row['predictionLog1pCPM'],predicted[j])
 for k,x in [('joint',predicted[j]),('noChange',c),('trainingMeanResponse',b)]:errors[k].append(x-y)
for k,e in errors.items():
 s=math.fsum(x*x for x in e);close(report['metrics'][k]['sumSquaredError'],s);close(report['metrics'][k]['RMSE'],math.sqrt(s/len(e)));close(report['metrics'][k]['MAE'],math.fsum(abs(x) for x in e)/len(e))
v={'status':'passed-independent-sparse-scoring','tag':tag,'compared':compared,'maximumScaledError':maximum,'scoreReportSHA256':hashlib.sha256((out/'scores.json').read_bytes()).hexdigest()};tmp=out/'scores-verification.tmp';tmp.write_text(json.dumps(v,indent=2)+'\n');tmp.replace(out/'scores-verification.json');print(v)
