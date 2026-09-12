from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import sys,json,hashlib,time
import numpy as np
from threadpoolctl import threadpool_limits
root=Path('/Users/n/numivivo-replogle2020-20260911');out=Path(__file__).parent
sys.path.insert(0,str(root))
from prepare_training import sha,dense_row,subset,encoded
protocol=json.loads((out/'protocol.json').read_text());grid=protocol['regularizationGrid']
assert grid==[0.01,0.1,1,10,100]
old=json.loads((root/'target-kernel/input-freeze.json').read_text())
annotations=json.loads((root/'identity-author-analysis/descriptors/annotations.json').read_text())
assert sha(root/'identity-author-analysis/descriptors/annotations.json')==old['annotationsSHA256']
contexts=json.loads((root/'training-preparation/count-folds.json').read_text())['contexts']
results=[];max_error=0.;t=time.time()
with threadpool_limits(limits=1):
 for context in contexts:
  gem=context['gemgroup'];path=root/'training-preparation'/(gem+'-training.json')
  assert sha(path)==context['trainingSHA256'];training=json.loads(path.read_text());targets=[x['id'] for x in training['selection']['targets']]
  for target in targets:
   keep=[i for i,x in enumerate(targets) if x!=target];selected=subset(training,keep)
   expected=next(x for x in old['folds'] if x['gemgroup']==gem and x['target']==target)
   selected_sha=hashlib.sha256(encoded(selected)).hexdigest();assert selected_sha==expected['selectedTrainingSHA256']
   # The selection function receives only the count-excluded training subset.
   names=[x['id'] for x in selected['selection']['targets']];assert len(names)==29 and target not in names
   counts=np.vstack([dense_row(selected['matrix'],i) for i in range(30)]).astype(float)
   expression=np.log1p(counts/counts.sum(axis=1,keepdims=True)*1e6);base=expression[0];y=expression[1:]-base
   terms=[set(annotations[x]['terms']) for x in names];k=np.array([[len(a&b)/len(a|b) for b in terms] for a in terms]);n=len(names)
   losses=[]
   for lam in grid:
    a=np.block([[k+lam*np.eye(n),np.ones((n,1))],[np.ones((1,n)),np.zeros((1,1))]])
    q=np.linalg.inv(a)[:n,:n];w=np.eye(n)-q/np.diag(q)[:,None]
    # Independent explicit leave-one-target-out solves verify all weights.
    explicit=np.zeros((n,n))
    for i in range(n):
     ix=[j for j in range(n) if j!=i];sub=k[np.ix_(ix,ix)]
     b=np.block([[sub+lam*np.eye(n-1),np.ones((n-1,1))],[np.ones((1,n-1)),np.zeros((1,1))]])
     explicit[i,ix]=np.linalg.solve(b,np.r_[k[i,ix],1])[:-1]
    error=float(np.max(abs(w-explicit)));max_error=max(max_error,error);assert error<1e-10
    prediction=np.maximum(base+w@y,0);rmse=np.sqrt(np.mean((prediction-expression[1:])**2,axis=1))
    losses.append(dict(regularization=lam,meanRMSE=float(rmse.mean()),innerTargets=names,innerRMSE=rmse.tolist()))
   # Exact equal losses choose the stronger regularization; no outer score used.
   chosen=min(losses,key=lambda x:(x['meanRMSE'],-x['regularization']))['regularization']
   results.append(dict(gemgroup=gem,target=target,selectedTrainingSHA256=selected_sha,selectedRegularization=chosen,losses=losses))
   (out/'progress.json').write_text(json.dumps(dict(folds=len(results),seconds=time.time()-t)))
assert len(results)==150
(out/'selection.json').write_text(json.dumps(dict(protocolSHA256=sha(out/'protocol.json'),selectorSHA256=sha(Path(__file__)),sourceInputFreezeSHA256=sha(root/'target-kernel/input-freeze.json'),maximumWeightOracleError=max_error,folds=results,outerScoringStarted=False,completedUnix=time.time()),indent=2))
print(json.dumps(dict(folds=len(results),seconds=time.time()-t,maximumWeightOracleError=max_error)))
