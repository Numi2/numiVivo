#!/usr/bin/env python3
"""Post-result matched fixed-lambda shuffle using immutable primary query weights."""
import argparse,json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from go_transfer import reference,select,independent_prediction
from combinations import counts,load,sha,write,metrics
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('mode',choices=['predict','score'])
for name in ['reference','primary','prepared','out']:p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--predictions',type=Path)
a=p.parse_args()
with threadpool_limits(limits=1):
 ref=reference(a.reference);pr=json.loads((a.primary/'receipt.json').read_text());dr=json.loads((a.prepared/'receipt.json').read_text())
 for path,key in [(a.primary/'predictions.npz','predictionsSHA256'),(a.primary/'folds.json','foldsSHA256')]:assert sha(path)==pr[key]
 assert sha(a.prepared/'kernel.npz')==dr['kernelSHA256']
 pred=load(a.primary/'predictions.npz');infos=json.loads((a.primary/'folds.json').read_text());desc=load(a.prepared/'kernel.npz');names=desc['targetIDs'].tolist();method='shuffledGoRidgeFixed'
 a.out.mkdir(parents=True,exist_ok=False)
 if a.mode=='predict':
  rows=[];errors=[]
  for i,target in enumerate(pred['targetIDs']):
   if not pred['supported'][i]:rows.append(np.full(len(pred['featureIDs']),np.nan));continue
   data=select(ref,target);training=infos[i]['descriptorTrainingTargets'];assert all(t in data['trainingTargets'] for t in training) and target not in training
   use=[data['trainingTargets'].index(t) for t in training];raw=counts(data['trainingCounts']);y=np.log1p(raw/raw.sum(axis=1,keepdims=True)*1e6)[use]-pred['baseline'][i];shuffled=np.roll(y,-1,axis=0)
   weight=np.array(infos[i]['weights']['goRidgeFixed']);prediction=weight@shuffled
   ix=[names.index(t) for t in training];kernel=desc['kernel'][np.ix_(ix,ix)];query=desc['kernel'][names.index(target),ix]
   oracle=independent_prediction(kernel,query,shuffled,1);np.testing.assert_allclose(prediction,oracle,atol=1e-10,rtol=1e-9);errors.append(float(np.max(abs(prediction-oracle))))
   rows.append(np.maximum(pred['baseline'][i]+prediction,0))
  np.savez_compressed(a.out/'predictions.npz',predictions=np.array(rows),targetIDs=pred['targetIDs'],featureIDs=pred['featureIDs'],supported=pred['supported'])
  write(a.out/'receipt.json',dict(predictionsSHA256=sha(a.out/'predictions.npz'),primaryReceiptSHA256=sha(a.primary/'receipt.json'),descriptorReceiptSHA256=sha(a.prepared/'receipt.json'),referenceSHA256=sha(a.reference),implementationSHA256=sha(Path(__file__)),protocolSHA256=sha(Path(__file__).with_name('GO_FIXED_CONTROL_PROTOCOL.md')),independentMaxAbsoluteError=max(errors),originalFixedWeightsUnchanged=True))
 else:
  if a.predictions is None:p.error('score requires --predictions')
  receipt=json.loads((a.predictions/'receipt.json').read_text());assert sha(a.predictions/'predictions.npz')==receipt['predictionsSHA256'];assert sha(a.primary/'receipt.json')==receipt['primaryReceiptSHA256'];control=load(a.predictions/'predictions.npz')
  np.testing.assert_array_equal(control['targetIDs'],pred['targetIDs']);np.testing.assert_array_equal(control['featureIDs'],pred['featureIDs']);np.testing.assert_array_equal(control['supported'],pred['supported'])
  results=[];conditions=ref['conditions'].tolist()
  for i,target in enumerate(pred['targetIDs']):
   if not pred['supported'][i]:continue
   raw=counts(ref['counts'][conditions.index(target)]);truth=np.log1p(raw/raw.sum()*1e6)-pred['baseline'][i]
   for panel,ix in [('allGenes',np.arange(len(truth))),('trainingTop1000',pred['trainingTop1000'][i])]:results.append(dict(target=str(target),method=method,panel=panel,**metrics((control['predictions'][i]-pred['baseline'][i])[ix],truth[ix])))
  write(a.out/'results.json',results);write(a.out/'summary.json',[dict(method=method,panel=panel,targets=sum(r['panel']==panel for r in results),meanResponseRMSE=float(np.mean([r['rmse'] for r in results if r['panel']==panel]))) for panel in ['allGenes','trainingTop1000']]);write(a.out/'receipt.json',dict(predictionReceiptSHA256=sha(a.predictions/'receipt.json'),referenceSHA256=sha(a.reference),implementationSHA256=sha(Path(__file__))))
