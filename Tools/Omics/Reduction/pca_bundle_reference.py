#!/usr/bin/env python3
"""Decode a verified PCA bundle for independent numerical comparison, not runtime I/O."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

def decode(root):
    model=json.loads((root/'model.json').read_text());metadata=json.loads((root/'metadata.json').read_text())
    quality=json.loads((root/'quality.json').read_text());receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['matrixFormat']=='complete-row-major-u32-row-u32-component-f64-le/v1'
    for key,name in [('model','model.json'),('metadata','metadata.json'),('quality','quality.json'),('scores','scores.bin'),('loadings','loadings.bin')]:
        with (root/name).open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==bytes(receipt[key]['bytes']).hex()
    rows=len(metadata['cells']);columns=model['options']['components'];selected=model['selectedFeatureIndices']
    assert rows==model['cells'] and len(quality)==rows
    def matrix(name, n):
        p=root/name;assert p.stat().st_size==n*columns*16
        records=np.fromfile(p,dtype=[('row','<u4'),('component','<u4'),('value','<f8')])
        np.testing.assert_array_equal(records['row'],np.repeat(np.arange(n),columns))
        np.testing.assert_array_equal(records['component'],np.tile(np.arange(columns),n))
        assert np.all(np.isfinite(records['value']))
        return records['value'].reshape(n,columns).copy()
    scores=matrix('scores.bin',rows);loadings=matrix('loadings.bin',len(selected))
    keys=['method','options','features','selectedFeatureIndices','explainedVariance','explainedVarianceRatio','relativeResiduals',
          'maximumLoadingOrthogonalityError','basisSize','projectionCenters','qualification']
    reduction={k:model[k] for k in keys};reduction['cells']=[dict(sampleID=c['sampleID'],barcode=c['barcode']) for c in metadata['cells']]
    reduction['scores']=scores.tolist();reduction['loadings']=loadings.tolist()
    return dict(reduction=reduction,reductionStorage=model['storage']),metadata,quality

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--baseline-report',type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    report,metadata,quality=decode(a.bundle)
    if a.baseline_report:
        old=json.loads(a.baseline_report.read_text());assert metadata==old['metadata'] and quality==old['quality']
        for key in ['scores','loadings','features','selectedFeatureIndices','explainedVariance','explainedVarianceRatio','relativeResiduals','maximumLoadingOrthogonalityError','basisSize','cells']:
            assert report['reduction'][key]==old['reduction'][key],key
        del old
    (a.out/'report.json').write_text(json.dumps(report,separators=(',',':'),allow_nan=False)+'\n')
    (a.out/'plan.json').write_bytes((a.bundle/'plan.json').read_bytes())
    (a.out/'bundle-checks.json').write_text(json.dumps(dict(status='passed',baselineFieldsExact=a.baseline_report is not None,cells=len(metadata['cells']),features=len(metadata['features']),allBinaryCoordinatesAndValuesChecked=True),indent=2)+'\n')
