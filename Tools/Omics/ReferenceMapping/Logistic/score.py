#!/usr/bin/env python3
"""Independent full-cohort optimizer and annotation assessment after prediction freeze."""
import argparse,importlib.metadata,json,tarfile,warnings
from pathlib import Path
import anndata as ad,numpy as np
from scipy.special import logsumexp,softmax
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score,f1_score,classification_report,confusion_matrix
from prepare import sha,write
from bundles import verify,read_json

def main():
 p=argparse.ArgumentParser();p.add_argument('--native',type=Path,required=True);p.add_argument('--prior',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 freeze=json.loads((a.native/'prediction-freeze.json').read_text());assert freeze['status']=='completed' and freeze['allFourDonors'] and len(freeze['archives'])==4
 rows=[]
 for item in freeze['archives']:
  donor=item['donor'];archive=a.native/item['path'];assert sha(archive)==item['SHA256'];members=verify(archive)
  with tarfile.open(archive,'r:gz') as t:
   model=read_json(t,members,'reference/model.json');report=read_json(t,members,'mapped/report.json')
  classes=np.array(model['classes']);labels=np.array(model['referenceLabels']);x=np.array(model['reduction']['scores']);q=np.array([c['scores'] for c in report['cells']]);native=np.array([c['classProbabilities'] for c in report['cells']]);pred=np.array([c['candidateLabel'] for c in report['cells']]);state=model['logistic']
  assert report['overlappingDonorIDs']==[] and report['distanceOperations']==0
  assert report['classifierOperations']==len(q)*(x.shape[1]+1)*len(classes)
  assert all(c['neighborIndices']==[] and c['squaredDistances']==[] and c['votes']==[] for c in report['cells'])
  scaler=StandardScaler().fit(x);tx=scaler.transform(x);qx=scaler.transform(q)
  np.testing.assert_allclose(state['means'],scaler.mean_,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(state['scales'],scaler.scale_,rtol=1e-12,atol=1e-12)
  y=np.searchsorted(classes,labels);counts=np.bincount(y,minlength=len(classes));weights=len(x)/(len(classes)*counts);np.testing.assert_array_equal(state['classCounts'],counts);np.testing.assert_allclose(state['classWeights'],weights,rtol=0,atol=0)
  params=np.array(state['parameters']);design=np.column_stack([tx,np.ones(len(tx))]);z=design@params.T;pt=softmax(z,axis=1);loss=np.sum(weights[y]*(logsumexp(z,axis=1)-z[np.arange(len(y)),y]))/len(y)+np.sum(params[:,:-1]**2)/(2*len(y))
  pt[np.arange(len(y)),y]-=1;gradient=(pt*(weights[y]/len(y))[:,None]).T@design;gradient[:,:-1]+=params[:,:-1]/len(y)
  assert abs(loss-state['objectives'][-1])<1e-11 and abs(np.abs(gradient).max()-state['gradientMaximum'])<1e-11
  assert np.abs(gradient).max()<=model['plan']['logistic']['gradientTolerance']+1e-11
  explicit=softmax(np.column_stack([qx,np.ones(len(qx))])@params.T,axis=1)
  np.testing.assert_allclose(native,explicit,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(native.sum(axis=1),1,rtol=0,atol=1e-14)
  np.testing.assert_array_equal(pred,classes[native.argmax(axis=1)])
  with warnings.catch_warnings(record=True) as caught:
   warnings.simplefilter('always');reference=LogisticRegression(C=1,solver='lbfgs',max_iter=10000,tol=1e-10,class_weight='balanced').fit(tx,labels)
  np.testing.assert_array_equal(reference.classes_,classes);rp=reference.predict_proba(qx);maximum=float(np.abs(native-rp).max());assert maximum<=1e-4,(donor,maximum)
  with np.load(a.prior/donor/'frozen-reference.npz',allow_pickle=False) as old:
   np.testing.assert_array_equal(model['featureIDs'],old['featureIDs']);np.testing.assert_array_equal(model['reduction']['selectedFeatureIndices'],old['selected']);np.testing.assert_array_equal(labels,old['referenceLabels'])
   np.testing.assert_array_equal([c['barcode'] for c in model['reduction']['cells']],old['referenceCells']);np.testing.assert_array_equal([c['cell']['barcode'] for c in report['cells']],old['queryCells'])
   trainError=float(np.abs(x-old['referenceScores']).max());queryError=float(np.abs(q-old['queryScores']).max());np.testing.assert_allclose(x,old['referenceScores'],rtol=1e-8,atol=1e-8);np.testing.assert_allclose(q,old['queryScores'],rtol=1e-8,atol=1e-8)
  source=ad.read_h5ad(a.prior/donor/'query/projected.h5ad',backed='r');truth=np.asarray(source.obs['cell_type'].astype(str));np.testing.assert_array_equal(source.obs_names,[c['cell']['barcode'] for c in report['cells']]);source.file.close()
  oldMetrics=json.loads((a.prior/donor/'metrics.json').read_text())['models']
  with np.load(a.prior/donor/'balancedLogistic-predictions.npz',allow_pickle=False) as old:
   np.testing.assert_array_equal(truth,old['truth']);priorProbabilityError=float(np.abs(native-old['probabilities']).max());priorDisagreements=int(np.sum(pred!=old['prediction']))
  per=classification_report(truth,pred,labels=classes,output_dict=True,zero_division=0);sortedp=np.sort(native,axis=1);disagreement=pred!=reference.predict(qx)
  result=dict(donor=donor,trainingCells=len(x),queryCells=len(q),classes=classes.tolist(),accuracy=float(accuracy_score(truth,pred)),macroF1=float(f1_score(truth,pred,labels=classes,average='macro',zero_division=0)),perLabel=per,confusion=confusion_matrix(truth,pred,labels=classes).tolist(),
   zeroRecallLabels=[c for c in classes if per[c]['support']>0 and per[c]['recall']==0],priorKNN={k:oldMetrics['knn15'][k] for k in ['accuracy','macroF1']},priorLogistic={k:oldMetrics['balancedLogistic'][k] for k in ['accuracy','macroF1']},
   independentProbabilityMaximumError=maximum,priorProbabilityMaximumError=priorProbabilityError,priorPredictionDisagreements=priorDisagreements,independentPredictionDisagreements=int(disagreement.sum()),independentDisagreementMargins=(sortedp[:,-1]-sortedp[:,-2])[disagreement].tolist(),independentWarnings=[str(w.message) for w in caught],
   independentGradientMaximum=float(np.abs(gradient).max()),independentObjective=float(loss),nativeIterations=state['iterations'],nativeEvaluations=state['evaluations'],chargedTrainingWork=state['chargedWork'],priorTrainingScoreMaximumError=trainError,priorQueryScoreMaximumError=queryError)
  rows.append(result);write(a.out/(donor+'.json'),result);np.savez_compressed(a.out/(donor+'-predictions.npz'),barcodes=np.asarray([c['cell']['barcode'] for c in report['cells']]),classes=classes,truth=truth,prediction=pred,probabilities=native,independentProbabilities=rp)
  print(json.dumps({k:result[k] for k in ['donor','queryCells','accuracy','macroF1','zeroRecallLabels','independentProbabilityMaximumError','priorPredictionDisagreements']}),flush=True)
 write(a.out/'summary.json',dict(status='passed',predictionFreezeSHA256=sha(a.native/'prediction-freeze.json'),checkerSHA256=sha(__file__),queryCells=sum(r['queryCells'] for r in rows),allDonors=4,folds=rows,packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','scikit-learn','anndata']},scope='Native balanced logistic annotation on inspected complete Baron donor folds; not calibrated probabilities, novel-class rejection, independent-study transfer or authoritative cell identities'))
if __name__=='__main__':main()
