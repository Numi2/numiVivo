#!/usr/bin/env python3
"""Check all complete Baron folds against Scanpy and the previously published native results."""
import argparse,json
from pathlib import Path
from pca_query_reference import check
from pca_bundle_reference import decode
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--source-folds',type=Path,required=True);p.add_argument('--prior',type=Path,required=True)
a=p.parse_args();results=[];identities=set()
for donor in ['human1','human2','human3','human4']:
    root=a.results/donor;source=a.source_folds/donor/'query/projected.h5ad';old=a.prior/donor
    checks,_=check(root/'query',source,root/'legacy-query-report.json')
    assert (root/'legacy-query-report.json').read_bytes()==(old/'mapped/report.json').read_bytes()
    current,training,_=decode(root/'query/reference');previous=json.loads((old/'reference/model.json').read_text())['reduction']
    for key in ['scores','loadings','features','selectedFeatureIndices','projectionCenters','explainedVariance','explainedVarianceRatio','relativeResiduals','maximumLoadingOrthogonalityError','basisSize','cells']:
        assert current['reduction'][key]==previous[key],(donor,key)
    metadata=json.loads((root/'query/metadata.json').read_text());report=json.loads((root/'query/report.json').read_text())
    assert report['overlappingDonorIDs']==[] and {s['donorID'] for s in metadata['samples']}=={donor}
    assert donor not in {s['donorID'] for s in training['samples']}
    cells={(c['sampleID'],c['barcode']) for c in metadata['cells']};assert not identities&cells;identities|=cells
    checks.update(donor=donor,legacyFullReportBytesExact=True,trainingFitNumericsExact=True,donorDisjoint=True)
    (root/'reference-checks.json').write_text(json.dumps(checks,indent=2)+'\n');results.append(checks);print(json.dumps(checks),flush=True)
assert len(identities)==8569
(a.results/'reference-checks.json').write_text(json.dumps(dict(status='passed',cells=len(identities),folds=results),indent=2)+'\n')
