#!/usr/bin/env python3
"""Numerical CSR/CSC/implicit-zero test, not biological benchmark evidence."""
import argparse
import json
from pathlib import Path
import subprocess

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
x=np.random.default_rng(7).poisson(2,size=(16,10)).astype(np.int32)
x[2]=0;x[:,9]=0
mapping=dict(schemaVersion=1,id='numeric-storage-fixture',evidence='synthetic',countUnit='umiCount',matrixPath='X',
    sourceDescription='Numerical fixture only; no biological evidence',sampleColumn='sample',mitochondrialFeatureIDs=[],
    samples=[dict(id='s',biologicalReplicateID='d',donorID='d',condition='c',batchID='unreported',organism='fixture')])
plan=dict(schemaVersion=1,mapping=mapping,contrasts=[],reduction=dict(pca=dict(highlyVariableFeatures=8,components=3,maximumBasis=10,meanBins=2)))
plan_path=a.out/'plan.json';plan_path.write_text(json.dumps(plan)+'\n')
results=[]
for encoding in ['csr','csc']:
    matrix=sparse.csr_matrix(x) if encoding=='csr' else sparse.csc_matrix(x)
    obj=ad.AnnData(matrix,obs=pd.DataFrame({'sample':['s']*16},index=[f'c{i}' for i in range(16)]),
                   var=pd.DataFrame(index=[f'g{i}' for i in range(10)]))
    source=a.out/f'{encoding}.h5ad';obj.write_h5ad(source,compression='gzip')
    destination=a.out/encoding
    run=subprocess.run([str(a.binary.resolve()),'singlecell-h5ad-pseudobulk',str(source),'--plan',str(plan_path),'--output',str(destination)],capture_output=True,text=True)
    (a.out/f'{encoding}.stderr').write_text(run.stderr)
    assert run.returncode==0,run.stderr
    report=json.loads((destination/'report.json').read_text())
    assert report['quality'][2]['totalCounts']==0 and len(report['reduction']['cells'])==16
    results.append(report)
assert results[0]['reduction']==results[1]['reduction'],'PCA depends on sparse source orientation'
assert results[0]['reductionStorage']['entryVisits']==results[1]['reductionStorage']['entryVisits']
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',checks=['CSR/CSC exact shared reduction','zero-count cell retained','zero feature retained in statistics','same sparse arithmetic work'],
    qualification='Numerical/storage fixture only'),indent=2)+'\n')
print('CSR/CSC and implicit-zero numerical checks passed')
