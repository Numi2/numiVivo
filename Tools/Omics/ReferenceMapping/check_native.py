#!/usr/bin/env python3
"""Qualify native frozen-reference mapping on the predeclared real Baron folds."""
import argparse, copy, hashlib, json, shutil, subprocess
from pathlib import Path
import anndata as ad
import numpy as np
from sklearn.neighbors import KNeighborsClassifier

p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--reference-results',type=Path,required=True);p.add_argument('--source-plan',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--scoped',action='store_true');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
plan=json.loads(a.source_plan.read_text());checks=[];commands=[]
def write(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def run(kind,args,log,fail=None):
    if a.scoped:
        command=[str(a.binary.resolve()),kind,*map(str,args)]
    elif kind in ['reference-verify','reference-map-verify']:
        command=[str(a.binary.resolve()),'singlecell-'+kind,*map(str,args)]
    elif kind=='reference-fit':
        source,plan,out=args;command=[str(a.binary.resolve()),'singlecell-reference-fit',str(source),'--plan',str(plan),'--output',str(out)]
    else:
        source,plan,ref,out=args;command=[str(a.binary.resolve()),'singlecell-reference-map',str(source),'--plan',str(plan),'--reference',str(ref),'--output',str(out)]
    r=subprocess.run(command,capture_output=True,text=True);log.write_text(r.stdout+r.stderr);commands.append({'command':command,'returncode':r.returncode,'expectedFailure':fail})
    write(a.out/'commands.json',commands)
    if fail:
        assert r.returncode!=0 and fail in r.stdout+r.stderr,(command,r.stdout,r.stderr)
    else:assert r.returncode==0,(command,r.stdout,r.stderr)
    checks.append(str(log.relative_to(a.out)))
summary=[]
for donor in ['human1','human2','human3','human4']:
    root=a.out/donor;root.mkdir();prior=a.reference_results/donor
    train=prior/'train/projected.h5ad';query=prior/'query/projected.h5ad'
    fit={'schemaVersion':1,'mapping':copy.deepcopy(plan['mapping']),'reduction':plan['reduction'],'featureNamespace':'HUMAN_GENE_SYMBOL','labelProvenance':'Baron 2016 source author cell_type assignments; no inferred labels; donor-held-out protocol.','neighbors':15}
    fit['mapping']['samples']=[s for s in fit['mapping']['samples'] if s['donorID']!=donor]
    qp={'schemaVersion':1,'mapping':copy.deepcopy(plan['mapping']),'featureNamespace':'HUMAN_GENE_SYMBOL'}
    qp['mapping'].pop('groupColumn',None);qp['mapping']['samples']=[s for s in qp['mapping']['samples'] if s['donorID']==donor]
    write(root/'fit.json',fit);write(root/'query.json',qp)
    run('reference-fit',[train,root/'fit.json',root/'reference'],root/'fit.log')
    run('reference-verify',[root/'reference'],root/'verify-reference.log')
    run('reference-map',[query,root/'query.json',root/'reference',root/'mapped'],root/'map.log')
    run('reference-map-verify',[root/'mapped'],root/'verify-map.log')
    run('reference-map',[query,root/'query.json',root/'reference',root/'repeat'],root/'repeat.log')
    assert (root/'mapped/report.json').read_bytes()==(root/'repeat/report.json').read_bytes()
    model=json.loads((root/'reference/model.json').read_text());report=json.loads((root/'mapped/report.json').read_text())
    with np.load(prior/'frozen-reference.npz',allow_pickle=False) as old:
        assert model['featureIDs']==old['featureIDs'].tolist()
        assert model['reduction']['selectedFeatureIndices']==old['selected'].tolist()
        assert model['referenceLabels']==old['referenceLabels'].tolist()
        assert [c['barcode'] for c in model['reduction']['cells']]==old['referenceCells'].tolist()
        np.testing.assert_allclose(model['reduction']['projectionCenters'],old['means'],rtol=1e-12,atol=1e-12)
        np.testing.assert_array_equal(model['reduction']['scores'],old['referenceScores'])
        scores=np.array([c['scores'] for c in report['cells']])
        assert [c['cell']['barcode'] for c in report['cells']]==old['queryCells'].tolist()
        np.testing.assert_allclose(scores,old['queryScores'],rtol=1e-9,atol=1e-10)
        maximumError=float(np.max(np.abs(scores-old['queryScores'])))
        reference=KNeighborsClassifier(n_neighbors=15,weights='uniform',algorithm='brute').fit(old['referenceScores'],old['referenceLabels'])
        distance,indices=reference.kneighbors(scores)
        np.testing.assert_array_equal([c['neighborIndices'] for c in report['cells']],indices)
        np.testing.assert_allclose([c['squaredDistances'] for c in report['cells']],distance**2,rtol=1e-9,atol=1e-9)
        np.testing.assert_array_equal([c['candidateLabel'] for c in report['cells']],reference.predict(scores))
        np.testing.assert_allclose(np.array([c['votes'] for c in report['cells']])/15,reference.predict_proba(scores),rtol=0,atol=0)
    assert report['overlappingDonorIDs']==[]
    with np.load(prior/'knn15-predictions.npz',allow_pickle=False) as old:
        np.testing.assert_array_equal([c['candidateLabel'] for c in report['cells']],old['prediction'])
    summary.append({'donor':donor,'cells':len(report['cells']),'maximumFrozenScoreError':maximumError,'exactNeighborIndicesAndLabels':True,'replay':True})
    print(json.dumps(summary[-1]),flush=True)
    if donor!='human1':continue
    # Mismatched contracts fail before mapping; existing outputs cannot be replaced.
    for name,change,error in [
        ('label-leak',lambda x:x['mapping'].update(groupColumn='cell_type'),'query label mapping is prohibited'),
        ('namespace',lambda x:x.update(featureNamespace='OTHER'),'namespace or count unit mismatch'),
        ('count-unit',lambda x:x['mapping'].update(countUnit='readCount'),'namespace or count unit mismatch'),
        ('distance-budget',lambda x:x.update(maximumDistanceOperations=1),'distance-operation budget'),
        ('projection-budget',lambda x:x.update(maximumProjectionUpdates=1),'projection-update budget'),
        ('organism',lambda x:[s.update(organism='NCBITaxon:10090') for s in x['mapping']['samples']],'gene universe and organism')]:
        bad=copy.deepcopy(qp);change(bad);write(root/(name+'.json'),bad)
        run('reference-map',[query,root/(name+'.json'),root/'reference',root/name],root/(name+'.log'),error)
        assert not (root/name).exists()
    run('reference-map',[query,root/'query.json',root/'reference',root/'mapped'],root/'overwrite.log','output already exists')
    # Query labels are removed completely and the gene axis reversed, while counts remain sparse.
    q=ad.read_h5ad(query);q.obs=q.obs.drop(columns=['cell_type']);q=q[:,::-1].copy();q.write_h5ad(root/'reordered-unlabelled.h5ad')
    run('reference-map',[root/'reordered-unlabelled.h5ad',root/'query.json',root/'reference',root/'reordered'],root/'reordered.log')
    other=json.loads((root/'reordered/report.json').read_text())
    np.testing.assert_allclose([c['scores'] for c in other['cells']],scores,rtol=1e-9,atol=1e-10)
    assert [c['candidateLabel'] for c in other['cells']]==[c['candidateLabel'] for c in report['cells']]
    q=q[:3].copy();q.X=q.X.tocsr();q.X.data[q.X.indptr[0]:q.X.indptr[1]]=0;q.X.eliminate_zeros();q.write_h5ad(root/'empty-library.h5ad')
    run('reference-map',[root/'empty-library.h5ad',root/'query.json',root/'reference',root/'empty'],root/'empty.log')
    empty=json.loads((root/'empty/report.json').read_text());assert empty['cells'][0].get('candidateLabel') is None and empty['cells'][0].get('scores') is None
    assert empty['distanceOperations']==2*len(model['referenceLabels'])*20
    q=q[:,1:].copy();q.write_h5ad(root/'missing-gene.h5ad')
    run('reference-map',[root/'missing-gene.h5ad',root/'query.json',root/'reference',root/'missing'],root/'missing.log','gene universe and organism')
    overlap=copy.deepcopy(qp);overlap['mapping']['samples']=fit['mapping']['samples'];write(root/'overlap.json',overlap)
    run('reference-map',[train,root/'overlap.json',root/'reference',root/'overlap'],root/'overlap.log','overlaps reference cell identities')
    shutil.copytree(root/'reference',root/'tampered-reference')
    path=root/'tampered-reference/model.json';path.write_bytes(path.read_bytes()+b' ')
    run('reference-verify',[root/'tampered-reference'],root/'tamper.log','bundle hash, schema or implementation mismatch')
    # A recomputed receipt hash must not turn a changed model into a valid fit.
    changed=json.loads(path.read_text());changed['referenceLabels'][0]='corrupted-label'
    write(path,changed)
    receiptPath=root/'tampered-reference/receipt.json';receipt=json.loads(receiptPath.read_text())
    receipt['result']={'bytes':list(hashlib.sha256(path.read_bytes()).digest())};write(receiptPath,receipt)
    run('reference-verify',[root/'tampered-reference'],root/'tamper-rehashed.log','reference model does not reconstruct')
    assert not list(root.glob('.numivivo-reference-*'))
write(a.out/'summary.json',{'folds':summary,'checks':checks,'binarySHA256':hashlib.file_digest(open(a.binary,'rb'),'sha256').hexdigest()})
