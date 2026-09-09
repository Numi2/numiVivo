#!/usr/bin/env python3
"""Check frozen query coordinates against independent Scanpy normalization/SciPy multiplication."""
import argparse,hashlib,json
from pathlib import Path
import anndata as ad
import numpy as np
import scanpy as sc
from pca_bundle_reference import decode

def check(root, source, prior=None):
    receipt=json.loads((root/'receipt.json').read_text());report=json.loads((root/'report.json').read_text())
    metadata=json.loads((root/'metadata.json').read_text());quality=json.loads((root/'quality.json').read_text())
    plan=json.loads((root/'plan.json').read_text());mapping=plan['mapping']
    for key,name in [('source',None),('plan','plan.json'),('reference','reference/receipt.json'),('metadata','metadata.json'),('quality','quality.json'),('report','report.json'),('scores','scores.bin')]:
        path=source if name is None else root/name
        with path.open('rb') as f: assert hashlib.file_digest(f,'sha256').digest()==bytes(receipt[key]['bytes']),name
    assert receipt['matrixFormat']=='complete-row-major-u32-row-u32-component-f64-le/v1'
    reduction,training,_=decode(root/'reference');r=reduction['reduction']
    trainplan=json.loads((root/'reference/plan.json').read_text())
    assert plan['featureNamespace']==trainplan['featureNamespace']
    n=len(metadata['cells']);d=r['options']['components']
    path=root/'scores.bin';assert path.stat().st_size==n*d*16
    records=np.fromfile(path,dtype=[('row','<u4'),('component','<u4'),('value','<f8')])
    np.testing.assert_array_equal(records['row'],np.repeat(np.arange(n),d));np.testing.assert_array_equal(records['component'],np.tile(np.arange(d),n))
    scores=records['value'].reshape(n,d);assert np.all(np.isfinite(scores))
    data=ad.read_h5ad(source)
    matrix=mapping['matrixPath'];x=data.X if matrix=='X' else data.layers[matrix.removeprefix('layers/')]
    x=x.tocsr().astype(np.float64);x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
    features=data.var_names.astype(str).tolist() if not mapping.get('featureIDColumn') else data.var[mapping['featureIDColumn']].astype(str).tolist()
    barcodes=data.obs_names.astype(str).tolist() if not mapping.get('barcodeColumn') else data.obs[mapping['barcodeColumn']].astype(str).tolist()
    assert features==[f['id'] for f in metadata['features']]
    assert barcodes==[c['barcode'] for c in metadata['cells']]
    assert data.obs[mapping['sampleColumn']].astype(str).tolist()==[c['sampleID'] for c in metadata['cells']]
    totals=np.asarray(x.sum(axis=1)).ravel()
    np.testing.assert_array_equal(totals,[q['totalCounts'] for q in quality]);np.testing.assert_array_equal(np.diff(x.indptr),[q['detectedFeatures'] for q in quality])
    observed=ad.AnnData(x);sc.pp.normalize_total(observed,target_sum=trainplan['reduction']['normalizationTarget']);sc.pp.log1p(observed)
    index={v:i for i,v in enumerate(features)};selected=[index[training['features'][j]['id']] for j in r['selectedFeatureIndices']]
    loadings=np.asarray(r['loadings']);centers=np.asarray(r['projectionCenters'])
    expected=observed.X[:,selected]@loadings-centers@loadings
    np.testing.assert_allclose(scores,expected,rtol=1e-9,atol=1e-10)
    assert report['selectedEntries']==observed.X[:,selected].nnz
    assert report['projectionUpdates']==report['selectedEntries']*d
    assert report['emptyLibraries']==int(np.sum(totals==0))
    assert not set((c['sampleID'],c['barcode']) for c in training['cells']) & set((c['sampleID'],c['barcode']) for c in metadata['cells'])
    checks=dict(status='passed',cells=n,features=len(features),components=d,maximumAbsoluteScoreError=float(np.max(np.abs(scores-expected))),allBinaryCoordinatesChecked=True,emptyLibraries=report['emptyLibraries'],projectionUpdates=report['projectionUpdates'])
    if prior:
        old=json.loads(prior.read_text())
        np.testing.assert_array_equal(scores,[c['scores'] for c in old['cells']])
        checks['legacyQueryScoresExact']=True
    return checks,scores

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--h5ad',type=Path,required=True);p.add_argument('--prior-report',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();checks,_=check(a.bundle,a.h5ad,a.prior_report);a.out.write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks))
