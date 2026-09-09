#!/usr/bin/env python3
"""Independent sparse Scanpy normalization and fixed-program numerical/biological checks."""
import argparse,hashlib,json
from pathlib import Path
import anndata as ad
import numpy as np
import scanpy as sc
from sklearn.metrics import roc_auc_score
p=argparse.ArgumentParser();p.add_argument('--prepared',type=Path,required=True);p.add_argument('--bundle',type=Path,required=True)
p.add_argument('--benchmark',choices=['kang','baron'],required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
report_path=a.bundle/'report.json';report=json.loads(report_path.read_text());result=report['programs']
plan=json.loads((a.bundle/'plan.json').read_text());mapping=plan['mapping']
source=a.prepared/'prepared.h5ad'
assert hashlib.sha256(source.read_bytes()).digest()==hashlib.sha256((a.bundle/'original.h5ad').read_bytes()).digest()
data=ad.read_h5ad(source)
assert mapping['matrixPath']=='X'
x=data.X.astype(np.float64).tocsr();x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
ids=data.var[mapping['featureIDColumn']].astype(str).tolist() if mapping.get('featureIDColumn') else data.var_names.astype(str).tolist()
assert len(set(ids))==len(ids)
barcodes=data.obs[mapping['barcodeColumn']].astype(str).tolist() if mapping.get('barcodeColumn') else data.obs_names.astype(str).tolist()
keys=list(zip(data.obs[mapping['sampleColumn']].astype(str),barcodes,strict=True));assert len(set(keys))==len(keys)
lookup={k:i for i,k in enumerate(keys)};order=np.array([lookup[(c['sampleID'],c['barcode'])] for c in result['cells']])
assert len(order)==len(keys) and len(set(order))==len(keys)
raw=x.copy();totals=np.asarray(x.sum(axis=1)).ravel()
norm=ad.AnnData(x);sc.pp.normalize_total(norm,target_sum=result['normalizationTarget']);sc.pp.log1p(norm)
x=norm.X.tocsr()[order];raw=raw[order]
reference=[];max_error=0.;programs=[]
for j,resolved in enumerate(result['programs']):
    definition=resolved['definition']
    assert definition==plan['programs']['definitions'][j]
    fingerprint=hashlib.sha256(json.dumps(definition,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).digest()
    assert fingerprint==bytes(resolved['definitionFingerprint']['bytes'])
    wanted={m['featureID']:m['weight'] for m in definition['members']}
    selected=[i for i,g in enumerate(ids) if g in wanted];missing=sorted(set(wanted)-set(ids))
    weights=np.array([wanted[ids[i]] for i in selected],dtype=np.float64)
    denominator=np.abs(weights).sum();weights/=denominator
    assert selected==resolved['featureIndices'] and missing==resolved['missingFeatureIDs']
    assert abs(denominator/sum(abs(w) for w in wanted.values())-resolved['weightCoverage'])<1e-14
    np.testing.assert_allclose(weights,resolved['effectiveWeights'],rtol=1e-14,atol=1e-15)
    expected=np.asarray(x[:,selected]@weights).ravel()
    observed=np.array([np.nan if row[j] is None else row[j] for row in result['scores']])
    empty=totals[order]==0;expected[empty]=np.nan
    np.testing.assert_allclose(observed,expected,rtol=2e-13,atol=2e-13,equal_nan=True)
    detected=raw[:,selected].getnnz(axis=1)
    np.testing.assert_array_equal(detected,np.array(result['detectedMembers'])[:,j])
    error=float(np.nanmax(np.abs(expected-observed)));max_error=max(max_error,error)
    reference.append(expected);programs.append(dict(id=definition['id'],requested=len(wanted),matched=len(selected),missing=missing,maximumAbsoluteScoreError=error))
# Predeclared descriptive observations: all eight paired donor IFN increases;
# beta-cell hallmark higher in author-labelled beta cells within each Baron donor.
samples={s['id']:s for s in report['metadata']['samples']}
donors=np.array([samples[c['sampleID']]['donorID'] for c in result['cells']])
conditions=np.array([samples[c['sampleID']]['condition'] for c in result['cells']])
biology=[]
if a.benchmark=='kang':
    assert len(set(donors))==8
    for j,values in enumerate(reference):
        for donor in sorted(set(donors)):
            ctrl=values[(donors==donor)&(conditions=='ctrl')];stim=values[(donors==donor)&(conditions=='stim')]
            assert len(ctrl)>0 and len(stim)>0
            delta=float(stim.mean()-ctrl.mean());biology.append(dict(program=programs[j]['id'],donor=donor,controlCells=len(ctrl),stimulatedCells=len(stim),meanDifference=delta))
    direction=all(r['meanDifference']>0 for r in biology)
else:
    assert len(set(donors))==4
    groups=np.array([str(data.obs.iloc[i][mapping['groupColumn']]) for i in order])
    for donor in sorted(set(donors)):
        mask=donors==donor;y=groups[mask]=='beta';v=reference[0][mask]
        assert y.any() and (~y).any()
        biology.append(dict(donor=donor,betaCells=int(y.sum()),otherCells=int((~y).sum()),meanDifference=float(v[y].mean()-v[~y].mean()),authorLabelAUC=float(roc_auc_score(y,v))))
    direction=all(r['meanDifference']>0 for r in biology)
out=dict(status='passed',numericalAgreement=True,expectedDirectionObserved=direction,cells=len(order),programs=programs,maximumAbsoluteScoreError=max_error,descriptiveBiology=biology,
    sourceSHA256=hashlib.sha256(source.read_bytes()).hexdigest(),reportSHA256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
    versions=dict(scanpy=sc.__version__,anndata=ad.__version__,numpy=np.__version__),
    qualification='Fixed expression-score agreement and within-donor descriptive observations. No automatic labels, inferred pathway activity, independent cell-type ground truth, held-out prediction or cross-study biological calibration.')
a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
