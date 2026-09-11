#!/usr/bin/env python3
"""Full sparse source-to-probability checks and explicitly coarse annotation diagnostics."""
import argparse,importlib.metadata,json,sys,tarfile,warnings
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd,scanpy as sc
from scipy.special import logsumexp,softmax
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score,balanced_accuracy_score,f1_score,classification_report,confusion_matrix
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'Logistic'))
from prepare import sha,write
from bundles import verify,read_json

FAMILIES={
 'kang':{'CD4 T cells':'T','CD8 T cells':'T','B cells':'B','NK cells':'NK','CD14+ Monocytes':'monocyte','FCGR3A+ Monocytes':'monocyte','Dendritic cells':'dendritic','Megakaryocytes':'megakaryocyte'},
 'ding':{'Cytotoxic T cell':'T','CD4+ T cell':'T','B cell':'B','Natural killer cell':'NK','CD14+ monocyte':'monocyte','CD16+ monocyte':'monocyte','Dendritic cell':'dendritic','Plasmacytoid dendritic cell':'dendritic','Megakaryocyte':'megakaryocyte'}}
FAMILY_ORDER=sorted(set(FAMILIES['kang'].values()))

def normal(a):
 x=a.X.astype(np.float64).tocsr();totals=np.asarray(x.sum(axis=1)).ravel();assert np.all(totals>0)
 x=x.multiply((10000/totals)[:,None]).tocsr();x.data=np.log1p(x.data)
 return x,totals

def metric(truth,pred,majority):
 return dict(cells=len(truth),accuracy=float(accuracy_score(truth,pred)),balancedAccuracy=float(balanced_accuracy_score(truth,pred)),macroF1=float(f1_score(truth,pred,labels=FAMILY_ORDER,average='macro',zero_division=0)),perFamily=classification_report(truth,pred,labels=FAMILY_ORDER,output_dict=True,zero_division=0),confusion=confusion_matrix(truth,pred,labels=FAMILY_ORDER).tolist(),trainingMajorityAccuracy=float(np.mean(truth==majority)))

def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--native',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 freeze=json.loads((a.native/'prediction-freeze.json').read_text());assert freeze['status']=='completed' and freeze['bothDirections'] and len(freeze['archives'])==2
 inputs=json.loads((a.inputs/'input-freeze.json').read_text())
 for name,h in inputs['files'].items():assert sha(a.inputs/name)==h,name
 results=[]
 for item in freeze['archives']:
  direction=item['direction'];train,query=direction.split('-to-');archive=a.native/item['path'];assert sha(archive)==item['SHA256'];members=verify(archive)
  with tarfile.open(archive,'r:gz') as t:
   model=read_json(t,members,'reference/model.json');report=read_json(t,members,'mapped/report.json')
  training=ad.read_h5ad(a.inputs/('kang.h5ad' if train=='kang' else 'ding-train.h5ad'))
  querying=ad.read_h5ad(a.inputs/('kang.h5ad' if query=='kang' else 'ding-query.h5ad'))
  tids=np.array(training.var_names if train=='kang' else training.var['reference_feature_id'],dtype=str)
  qids=np.array(querying.var_names if query=='kang' else querying.var['reference_feature_id'],dtype=str)
  np.testing.assert_array_equal(tids,model['featureIDs']);np.testing.assert_array_equal(training.obs_names,[c['barcode'] for c in model['reduction']['cells']]);np.testing.assert_array_equal(querying.obs_names,[c['cell']['barcode'] for c in report['cells']])
  np.testing.assert_array_equal(training.obs['native_sample'].astype(str),[c['sampleID'] for c in model['reduction']['cells']]);np.testing.assert_array_equal(querying.obs['native_sample'].astype(str),[c['cell']['sampleID'] for c in report['cells']])
  labels=np.asarray(training.obs['cell_type' if train=='kang' else 'sourceCellType'].astype(str));np.testing.assert_array_equal(labels,model['referenceLabels']);classes=np.array(model['classes']);assert set(classes)==set(FAMILIES[train])
  x,totals=normal(training);q,qtotals=normal(querying);np.testing.assert_array_equal(qtotals,[c['totalCounts'] for c in report['cells']])
  reduction=model['reduction'];panel=set(model['plan']['reduction']['pca']['featurePanel']);pi=np.array([i for i,g in enumerate(tids) if g in panel]);assert len(pi)==len(panel)==14976
  selected=np.array(reduction['selectedFeatureIndices']);assert set(tids[selected]).issubset(panel)
  hvg=ad.AnnData(x[:,pi].copy());hvg.var_names=tids[pi]
  sc.pp.highly_variable_genes(hvg,flavor='seurat',n_top_genes=2000,n_bins=20)
  np.testing.assert_array_equal(selected,pi[np.asarray(hvg.var['highly_variable'])])
  means=np.asarray(x[:,selected].mean(axis=0)).ravel();v=np.array(reduction['loadings']);nativeX=np.array(reduction['scores']);nativeQ=np.array([c['scores'] for c in report['cells']]);variance=np.array(reduction['explainedVariance'])
  np.testing.assert_allclose(means,reduction['projectionCenters'],rtol=1e-11,atol=1e-11)
  sx=x[:,selected];directX=np.asarray(sx@v)-means@v;np.testing.assert_allclose(directX,nativeX,rtol=1e-7,atol=1e-7)
  covarianceV=(sx.T@directX-means[:,None]*directX.sum(axis=0)[None,:])/(len(training)-1)
  residual=np.linalg.norm(covarianceV-v*variance[None,:],axis=0)/variance;assert residual.max()<=1.01e-6,float(residual.max())
  qi={g:i for i,g in enumerate(qids)};directQ=np.asarray(q[:,[qi[g] for g in tids[selected]]]@v)-means@v;np.testing.assert_allclose(directQ,nativeQ,rtol=1e-7,atol=1e-7)
  assert report['distanceOperations']==0 and report['classifierOperations']==len(querying)*21*len(classes)
  assert all(c['neighborIndices']==[] and c['squaredDistances']==[] and c['votes']==[] for c in report['cells'])
  state=model['logistic'];scaler=StandardScaler().fit(nativeX);tx=scaler.transform(nativeX);qx=scaler.transform(nativeQ)
  np.testing.assert_allclose(state['means'],scaler.mean_,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(state['scales'],scaler.scale_,rtol=1e-12,atol=1e-12)
  y=np.searchsorted(classes,labels);counts=np.bincount(y,minlength=len(classes));weights=len(tx)/(len(classes)*counts);np.testing.assert_array_equal(state['classCounts'],counts);np.testing.assert_allclose(state['classWeights'],weights,rtol=0,atol=0)
  params=np.array(state['parameters']);design=np.column_stack([tx,np.ones(len(tx))]);z=design@params.T;pt=softmax(z,axis=1)
  objective=float(np.sum(weights[y]*(logsumexp(z,axis=1)-z[np.arange(len(y)),y]))/len(y)+np.sum(params[:,:-1]**2)/(2*len(y)))
  pt[np.arange(len(y)),y]-=1;gradient=(pt*(weights[y]/len(y))[:,None]).T@design;gradient[:,:-1]+=params[:,:-1]/len(y)
  assert abs(objective-state['objectives'][-1])<1e-11;assert np.abs(gradient).max()<=model['plan']['logistic']['gradientTolerance']+1e-11
  nativeP=np.array([c['classProbabilities'] for c in report['cells']]);explicit=softmax(np.column_stack([qx,np.ones(len(qx))])@params.T,axis=1);np.testing.assert_allclose(nativeP,explicit,rtol=1e-12,atol=1e-12)
  pred=np.array([c['candidateLabel'] for c in report['cells']]);np.testing.assert_array_equal(pred,classes[nativeP.argmax(axis=1)])
  with warnings.catch_warnings(record=True) as caught:
   warnings.simplefilter('always');reference=LogisticRegression(C=1,solver='lbfgs',class_weight='balanced',max_iter=10000,tol=1e-10).fit(tx,labels)
  rp=reference.predict_proba(qx);maximum=float(np.abs(nativeP-rp).max());assert maximum<=1e-4,(direction,maximum)
  # Source query labels enter only after every source/numerical comparison passes.
  truth=np.asarray(querying.obs['cell_type' if query=='kang' else 'sourceCellType'].astype(str));assigned=np.array([t in FAMILIES[query] for t in truth]);tf=np.array([FAMILIES[query].get(t,'unavailable') for t in truth]);pf=np.array([FAMILIES[train][t] for t in pred]);majority=FAMILIES[train][classes[counts.argmax()]]
  metrics=metric(tf[assigned],pf[assigned],majority);byLibrary={}
  for library in sorted(set(querying.obs['native_sample'].astype(str))):
   mask=np.asarray(querying.obs['native_sample'].astype(str)==library);valid=mask&assigned
   byLibrary[library]=dict(queryCells=int(mask.sum()),assignedCells=int(valid.sum()),metrics=metric(tf[valid],pf[valid],majority) if valid.any() else None)
  table=pd.crosstab(pd.Series(truth,name='sourceLabel'),pd.Series(pred,name='nativeCandidateLabel')).reindex(columns=classes,fill_value=0)
  minimumRecall=min(metrics['perFamily'][f]['recall'] for f in FAMILY_ORDER if metrics['perFamily'][f]['support']>0)
  result=dict(direction=direction,trainingCells=len(training),queryCells=len(querying),assignedQueryCells=int(assigned.sum()),unassignedQueryCells=int((truth=='Unassigned').sum()),unavailableQueryCells=int((truth=='').sum()),panelFeatures=len(panel),selectedFeatures=len(selected),sourceNormalization='Complete measured RNA gene universe separately for each source; shared-panel subtotal is never substituted',HVGSelectionExact=True,maximumPCACovarianceResidual=float(residual.max()),maximumTrainingProjectionError=float(np.abs(directX-nativeX).max()),maximumQueryProjectionError=float(np.abs(directQ-nativeQ).max()),maximumIndependentProbabilityError=maximum,independentLabelDisagreements=int(np.sum(pred!=reference.predict(qx))),independentGradientMaximum=float(np.abs(gradient).max()),independentObjective=objective,independentWarnings=[str(w.message) for w in caught],originalPredictionClasses=classes.tolist(),originalSourceLabels=table.index.tolist(),originalRectangularConfusion=table.to_numpy().tolist(),familyOrder=FAMILY_ORDER,coarseFamilyMetrics=metrics,coarseFamilyTargetPassed=bool(metrics['balancedAccuracy']>=.8 and metrics['macroF1']>=.8 and minimumRecall>=.5),byLibrary=byLibrary,scope='Complete reused-cohort transfer and coarse family diagnostic; original fine labels are not a shared taxonomy, source assignments are not experimental identity proof, probabilities uncalibrated.')
  write(a.out/(direction+'.json'),result);results.append(result)
  np.savez_compressed(a.out/(direction+'-predictions.npz'),barcodes=np.asarray(querying.obs_names,dtype=str),sampleIDs=np.asarray(querying.obs['native_sample'],dtype=str),classes=classes,sourceLabels=truth,predictions=pred,probabilities=nativeP,independentProbabilities=rp,sourceFamilies=tf,predictedFamilies=pf)
  print(json.dumps({k:result[k] for k in ['direction','queryCells','assignedQueryCells','maximumIndependentProbabilityError','coarseFamilyTargetPassed']}),flush=True)
 write(a.out/'summary.json',dict(status='passed-numerical-checks',predictionFreezeSHA256=sha(a.native/'prediction-freeze.json'),checkerSHA256=sha(__file__),packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','scikit-learn','anndata','scanpy']},queryCells=sum(r['queryCells'] for r in results),directions=results))
if __name__=='__main__':main()
