#!/usr/bin/env python3
"""Adapt verified file-backed graphs to the independent exact-kNN/umap-learn checker."""
import argparse,gzip,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
from pca_bundle_reference import decode
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--baseline-pca-report',type=Path)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
receipt=json.loads((a.bundle/'receipt.json').read_text());graph=json.loads((a.bundle/'graph.json').read_text());plan=json.loads((a.bundle/'plan.json').read_text())
for key,name in [('input','input/receipt.json'),('plan','plan.json'),('graph','graph.json'),('executionReport','execution.json')]:
    with (a.bundle/name).open('rb') as f:assert hashlib.file_digest(f,'sha256').digest()==bytes(receipt[key]['bytes']),name
source=a.bundle/'input';ir=json.loads((source/'receipt.json').read_text());metadata=json.loads((source/'metadata.json').read_text())
for key,name in [('scores','scores.bin'),('metadata','metadata.json')]:
    with (source/name).open('rb') as f:assert hashlib.file_digest(f,'sha256').digest()==bytes(ir[key]['bytes']),name
cells=[dict(sampleID=c['sampleID'],barcode=c['barcode']) for c in metadata['cells']];assert cells==graph['cells']
n=len(cells);d=graph['dimensions'];path=source/'scores.bin';assert path.stat().st_size==n*d*16
records=np.fromfile(path,dtype=[('row','<u4'),('component','<u4'),('value','<f8')]);np.testing.assert_array_equal(records['row'],np.repeat(np.arange(n),d));np.testing.assert_array_equal(records['component'],np.tile(np.arange(d),n));assert np.all(np.isfinite(records['value']))
execution=json.loads((a.bundle/'execution.json').read_text());assert execution['directedDistanceEvaluations']==n*(n-1) and execution['scalarDistanceTerms']==n*(n-1)*d
assert execution['scoreRecordReads']==n*d*((n+plan['execution']['queryBlockRows']-1)//plan['execution']['queryBlockRows']+1)
if a.baseline_pca_report:
    current,_,_=decode(source)
    data=a.baseline_pca_report.read_bytes();old=json.loads(gzip.decompress(data) if a.baseline_pca_report.suffix=='.gz' else data)['reduction']
    for key in ['scores','loadings','features','selectedFeatureIndices','explainedVariance','explainedVarianceRatio','relativeResiduals','maximumLoadingOrthogonalityError','basisSize','cells']:assert current['reduction'][key]==old[key],key
report=dict(neighbors=graph,reduction=dict(cells=cells,scores=records['value'].reshape(n,d).tolist()))
(a.out/'report.json').write_text(json.dumps(report,separators=(',',':'),allow_nan=False)+'\n')
r=subprocess.run([sys.executable,str(Path(__file__).with_name('check_neighbors.py')),'--report',str(a.out/'report.json'),'--out',str(a.out/'checks.json')],capture_output=True,text=True)
(a.out/'run.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stdout+r.stderr
checks=json.loads((a.out/'checks.json').read_text());checks.update(allBinaryCoordinatesChecked=True,executionAccountingChecked=True,baselinePCANumericsExact=a.baseline_pca_report is not None)
(a.out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks))
