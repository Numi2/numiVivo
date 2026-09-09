#!/usr/bin/env python3
"""Haber 2017: all annotated tuft cells, native counts/QC versus Scanpy.

No DE claim: treatment groups have only two batch labels and individual donor
identity is not supplied. Batch labels are not invented biological replicates.
The source dense matrix is read in 64-row blocks into sparse storage.
"""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import time
import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
source_sha = 'efcb49770529ec14e1ba6b45b3d0086b3829e4a63a85c02793f0d7ad90b6db21'
assert hashlib.sha256(a.source.read_bytes()).hexdigest() == source_sha
source = ad.read_h5ad(a.source, backed='r')
assert source.shape == (9842,15215)
rows = np.flatnonzero((source.obs.cell_label == 'Tuft').to_numpy())
counts = sparse.vstack([sparse.csr_matrix(source[rows[i:i+64],:].X) for i in range(0,len(rows),64)],format='csr')
assert counts.shape == (409,15215) and counts.nnz == 684794
assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and np.equal(counts.data,np.floor(counts.data)).all()
counts = counts.astype(np.int64)
obs,var = source.obs.iloc[rows].copy(),source.var.copy()
source.file.close()
assert obs.index.is_unique and var.index.is_unique
prepared = ad.AnnData(X=counts,obs=obs,var=var)
prepared.uns['numivivo_scope'] = dict(source_sha256=source_sha, source_rows=rows, scope='All annotated tuft cells and all source genes; QC only')
prepared.write_h5ad(a.out/'prepared.h5ad',compression='gzip')
samples = []
for batch,group in obs.groupby('batch',observed=True):
    assert group.condition.nunique() == 1
    samples.append(dict(id=str(batch), biologicalReplicateID='unreported', batchID=str(batch),
        condition=str(group.condition.iloc[0]), organism='NCBITaxon:10090'))
def save(name,obj):
    (a.out/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
save('mapping.json',dict(schemaVersion=1,id='haber2017-tuft-qc',evidence='measured',countUnit='umiCount',matrixPath='X',
    sampleColumn='batch',groupColumn='cell_label',mitochondrialFeatureIDs=[],samples=samples,
    sourceDescription='Haber 2017 doi:10.1038/nature24489, all annotated tuft cells; no donor-level inference'))
save('analysis.json',dict(schemaVersion=1,id='haber-tuft-qc',normalizationTarget=10000,contrasts=[]))
commands=[]
def run(label,*args):
    cmd=[str(a.binary.resolve()),*map(str,args)]
    start=time.perf_counter()
    result=subprocess.run(['/usr/bin/time','-l',*cmd],capture_output=True,text=True)
    (a.out/(label+'.log')).write_text(result.stdout+result.stderr)
    peak=next((int(line.split()[0]) for line in result.stderr.splitlines() if 'maximum resident set size' in line),None)
    commands.append(dict(label=label,command=cmd,returncode=result.returncode,elapsedSeconds=time.perf_counter()-start,peakResidentBytes=peak))
    save('commands.json',commands)
    assert result.returncode==0,result.stderr
run('import','singlecell-h5ad-import',a.out/'prepared.h5ad','--plan',a.out/'mapping.json','--output',a.out/'imported')
run('counts','singlecell-run',a.out/'imported/manifest.json','--store',a.out/'store','--output',a.out/'counts.json')
run('analysis','singlecell-analyze',a.out/'counts.json','--plan',a.out/'analysis.json','--store',a.out/'store','--output',a.out/'analysis-receipt.json')
run('replay','singlecell-analysis-export',a.out/'analysis-receipt.json','--store',a.out/'store','--output',a.out/'native-report.json')
native=json.loads((a.out/'native-report.json').read_text())['processed']
index={str(v):i for i,v in enumerate(obs.index)}
order=np.array([index[c['barcode']] for c in native['dataset']['cells']])
assert len(order)==len(rows) and len(set(order))==len(rows)
m=native['dataset']['matrix']
recovered=sparse.csr_matrix((m['counts'],m['featureIndices'],m['rowOffsets']),shape=counts.shape)
assert (recovered!=counts[order]).nnz==0
reference=ad.AnnData(X=counts.astype(np.float64).copy())
sc.pp.calculate_qc_metrics(reference,percent_top=None,log1p=False,inplace=True)
quality=[d['quality'] for d in native['decisions']]
np.testing.assert_array_equal([q['totalCounts'] for q in quality],reference.obs.total_counts.to_numpy()[order])
np.testing.assert_array_equal([q['detectedFeatures'] for q in quality],reference.obs.n_genes_by_counts.to_numpy()[order])
sc.pp.normalize_total(reference,target_sum=10000)
sc.pp.log1p(reference)
m=native['normalized']
normalized=sparse.csr_matrix((m['values'],m['featureIndices'],m['rowOffsets']),shape=counts.shape)
r=reference.X[order].tocsr()
np.testing.assert_array_equal(normalized.indices,r.indices)
np.testing.assert_array_equal(normalized.indptr,r.indptr)
error=float(np.max(np.abs(normalized.data-r.data)))
assert error<1e-10
report=dict(schemaVersion=1,status='passed-real-counts-qc',sourceSHA256=source_sha,
    sourceURL='https://ndownloader.figshare.com/files/54169301',publication='https://doi.org/10.1038/nature24489',
    scope='All annotated tuft cells and all source genes',cells=len(rows),features=counts.shape[1],nonzeros=counts.nnz,
    exactImportedCounts=True,exactScanpyQC=True,maxLogNormalizationError=error,
    donorDEQualified=False,limitation='Individual donors are not identified; two batch labels per treatment are insufficient for the current three-replicate gate',
    versions={k:version(k) for k in ['scanpy','anndata','numpy','scipy','pandas','h5py']},
    binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),nativeCommands=commands)
save('report.json',report)
print(json.dumps({k:v for k,v in report.items() if k!='nativeCommands'},indent=2))
