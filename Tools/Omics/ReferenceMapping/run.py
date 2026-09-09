#!/usr/bin/env python3
"""Run the predeclared Baron donor-held-out reference benchmark; sparse counts only."""
import argparse, hashlib, importlib.metadata, json, subprocess, warnings
from pathlib import Path
import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, confusion_matrix, f1_score


def digest(p):
    return hashlib.file_digest(open(p, 'rb'), 'sha256').hexdigest()


def write(p, value):
    p.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def project(x, selected, means, loadings):
    # This function deliberately has no query labels or fitting operation.
    return np.asarray(x[:, selected] @ loadings) - means @ loadings


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=False)
    source_hash=digest(a.source)
    write(a.out/'inputs.json',{'sourceSHA256':source_hash,'binarySHA256':digest(a.binary),
          'planSHA256':digest(a.plan),'packages':{n:importlib.metadata.version(n) for n in ['anndata','scanpy','numpy','scipy','scikit-learn']}})
    obj=ad.read_h5ad(a.source)
    assert sparse.issparse(obj.X) and obj.var_names.is_unique and obj.obs_names.is_unique
    donors=obj.obs['donor'].astype(str).to_numpy()
    labels=obj.obs['cell_type'].astype(str).to_numpy()
    classes=sorted(set(labels)); results=[]
    def run(args,log):
        with log.open('w') as f:
            subprocess.run([str(a.binary.resolve()),*map(str,args)],stdout=f,stderr=subprocess.STDOUT,check=True)
    for donor in sorted(set(donors)):
        fold=a.out/donor; fold.mkdir()
        train=np.flatnonzero(donors!=donor); query=np.flatnonzero(donors==donor)
        for name,rows in [('train',train),('query',query)]:
            plan={'schemaVersion':1,'source':{'bytes':list(bytes.fromhex(source_hash))},
                  'provenance':'Baron donor-held-out reference mapping protocol; '+donor+' '+name,
                  'observationIndices':rows.tolist()}
            write(fold/(name+'.json'),plan)
            run(['singlecell-h5ad-project',a.source,'--plan',fold/(name+'.json'),'--output',fold/name],fold/(name+'.log'))
        plan=json.loads(a.plan.read_text()); plan.pop('programs',None)
        plan['mapping']['samples']=[s for s in plan['mapping']['samples'] if s['donorID']!=donor]
        write(fold/'analysis.json',plan)
        run(['singlecell-h5ad-pseudobulk',fold/'train/projected.h5ad','--plan',fold/'analysis.json','--output',fold/'native'],fold/'native.log')
        report=json.loads((fold/'native/report.json').read_text()); r=report['reduction']
        b=ad.read_h5ad(fold/'train/projected.h5ad'); q=ad.read_h5ad(fold/'query/projected.h5ad')
        assert list(b.var_names)==list(q.var_names)==list(obj.var_names)
        # Match native reference row order explicitly; source author labels never enter PCA.
        mapping=plan['mapping']; sample=mapping['sampleColumn']
        lookup={(str(s),str(c)):i for i,(s,c) in enumerate(zip(b.obs[sample],b.obs_names))}
        order=[lookup[(c['sampleID'],c['barcode'])] for c in r['cells']]
        assert len(order)==len(b) and len(set(order))==len(b)
        assert set(b.obs_names).isdisjoint(q.obs_names)
        assert set(b.obs['donor'].astype(str))==set(donors)-{donor}
        assert set(q.obs['donor'].astype(str))=={donor}
        b=b[order].copy(); y=np.asarray(b.obs['cell_type'].astype(str),dtype=str)
        for data in [b,q]:
            data.X=data.X.astype(np.float64).tocsr()
            sc.pp.normalize_total(data,target_sum=10000); sc.pp.log1p(data)
            assert sparse.issparse(data.X)
        sc.pp.highly_variable_genes(b,flavor='seurat',n_top_genes=2000,n_bins=20)
        selected=np.flatnonzero(b.var['highly_variable'].to_numpy())
        assert selected.tolist()==r['selectedFeatureIndices']
        ref=b[:,selected].copy()
        sc.pp.pca(ref,n_comps=20,zero_center=True,svd_solver='arpack',dtype='float64',random_state=7)
        np.testing.assert_allclose(r['explainedVariance'],ref.uns['pca']['variance'],rtol=1e-7,atol=1e-9)
        means=np.asarray(b.X[:,selected].mean(axis=0)).ravel()
        loadings=np.asarray(r['loadings']); scores=np.asarray(r['scores'])
        np.testing.assert_allclose(project(b.X,selected,means,loadings),scores,rtol=1e-7,atol=1e-8)
        projected=project(q.X,selected,means,loadings)
        scaler=StandardScaler().fit(scores)
        np.savez_compressed(fold/'frozen-reference.npz',featureIDs=np.asarray(b.var_names,dtype=str),selected=selected,means=means,loadings=loadings,
                            referenceScores=scores,referenceLabels=y,referenceCells=np.asarray(b.obs_names,dtype=str),queryCells=np.asarray(q.obs_names,dtype=str),
                            queryScores=projected,scaleMean=scaler.mean_,scale=scaler.scale_)
        metrics={'heldOutDonor':donor,'trainingCells':len(train),'queryCells':len(query),'absentTrainingLabels':sorted(set(classes)-set(y)),'models':{}}
        for name,model,tx,qx in [
            ('knn15',KNeighborsClassifier(n_neighbors=15,weights='uniform',algorithm='brute'),scores,projected),
            ('balancedLogistic',LogisticRegression(C=1,solver='lbfgs',max_iter=2000,tol=1e-6,class_weight='balanced'),scaler.transform(scores),scaler.transform(projected))]:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always'); model.fit(tx,y); pred=model.predict(qx); probs=model.predict_proba(qx)
            truth=np.asarray(q.obs['cell_type'].astype(str),dtype=str)
            metrics['models'][name]={'accuracy':accuracy_score(truth,pred),'balancedAccuracy':balanced_accuracy_score(truth,pred),
                'macroF1':f1_score(truth,pred,labels=classes,average='macro',zero_division=0),
                'perLabel':classification_report(truth,pred,labels=classes,output_dict=True,zero_division=0),
                'classOrder':classes,'confusion':confusion_matrix(truth,pred,labels=classes).tolist(),
                'warnings':[str(w.message) for w in caught],'parameters':model.get_params()}
            np.savez_compressed(fold/(name+'-predictions.npz'),prediction=pred,truth=truth,probabilities=probs,classes=model.classes_)
            if hasattr(model,'coef_'):
                np.savez_compressed(fold/(name+'-model.npz'),coefficients=model.coef_,intercept=model.intercept_,classes=model.classes_,iterations=model.n_iter_)
        write(fold/'metrics.json',metrics); results.append(metrics)
        print(json.dumps({'donor':donor,'scores':{n:{k:v[k] for k in ['accuracy','balancedAccuracy','macroF1']} for n,v in metrics['models'].items()}}),flush=True)
    write(a.out/'results.json',results)

if __name__=='__main__': main()
