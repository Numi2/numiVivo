#!/usr/bin/env python3
"""Prepare all source Kang cells; reuse the published, source-supported Hagai design."""
import argparse,hashlib,json
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse
p=argparse.ArgumentParser(description=__doc__)
for name in ['source','hagai-fit','baron-fit','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
repo=Path(__file__).resolve().parents[3]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(name,v):
 p=a.out/name;p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');return p
assert sha(a.source)=='e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830'
x=ad.read_h5ad(a.source);assert x.shape==(24673,15706)
counts=sparse.csr_matrix(x.X);assert np.isfinite(counts.data).all() and (counts.data>=0).all() and np.equal(counts.data,np.floor(counts.data)).all()
obs=x.obs[['cell_type','replicate','label']].copy();obs['native_sample']=obs.replicate.astype(str)+'__'+obs.label.astype(str)
y=ad.AnnData(X=counts.astype(np.int64),obs=obs,var=x.var.copy());y.write_h5ad(a.out/'kang.h5ad',compression='gzip')
loaded=ad.read_h5ad(a.out/'kang.h5ad');assert (loaded.X!=counts).nnz==0 and loaded.obs_names.equals(x.obs_names) and loaded.var_names.equals(x.var_names)
pairs=sorted(set(zip(obs.replicate.astype(str),obs.label.astype(str))));assert len(pairs)==16
mapping=json.loads((repo/'Tools/Omics/Benchmarks/evidence/2026-09-09-hagai-nb/plan.json').read_text())['mapping']
base=json.loads(a.hagai_fit.read_text());base['mapping']=mapping;save('hagai-fit.json',base)
base['mapping']=dict(schemaVersion=1,id='kang2018-all-cells',evidence='measured',countUnit='umiCount',matrixPath='X',
 sourceDescription='Kang 2018 doi:10.1038/nbt.4042, pertpy file 34464122; all 24673 cells and 15706 genes in source order. Eight paired donors, eight source-annotated cell types. Batch unreported.',
 sampleColumn='native_sample',groupColumn='cell_type',mitochondrialFeatureIDs=[],samples=[dict(id=d+'__'+c,biologicalReplicateID=d,donorID=d,condition=c,batchID='unreported',organism='NCBITaxon:9606') for d,c in pairs])
save('kang-fit.json',base);save('baron-fit.json',json.loads(a.baron_fit.read_text()))
composition=obs.groupby(['replicate','cell_type','label'],observed=True).size().reset_index(name='cells').astype({'replicate':str,'cell_type':str,'label':str}).to_dict('records')
save('preparation.json',dict(sourceURL='https://ndownloader.figshare.com/files/34464122',sourceSHA256=sha(a.source),preparedSHA256=sha(a.out/'kang.h5ad'),cells=x.n_obs,features=x.n_vars,nonzeros=counts.nnz,allCountValuesExact=True,allCellFeatureIdentitiesExact=True,composition=composition,hagaiDesignSource='Tools/Omics/Benchmarks/evidence/2026-09-09-hagai-nb/plan.json',hagaiPairing='Prefixes mouse1/2/3 identify paired timecourses; Supplementary Table 2 supports three individuals, no prefix-to-table-row mapping claimed'))
protocol=dict(schemaVersion=1,cohorts=dict(kang=24673,hagai=13863,baron=8569),seeds=[7,19,41],integration=dict(covariate='donor',clusters=100,diversity=2,ridge=1,temperature=0.1,maximumIterations=10,relativeTolerance=0.01,maximumWork=1000000000),evaluationNeighbors=30,
 gateMargins=dict(maximumCellTypeBalancedAccuracyLoss=0.02,maximumPerTypeRecallLoss=0.05,maximumConditionBalancedAccuracyLoss=0.02,maximumProgramSpearmanLoss=0.05,minimumNegativeControlAccuracyLoss=0.10),
 program=dict(kang=['IFI6','IFIT1','ISG15','MX1','ISG20'],hagai=['ENSMUSG00000025225','ENSMUSG00000021025','ENSMUSG00000034855','ENSMUSG00000035692','ENSMUSG00000024401']),
 scope='Full cohorts, all genes at import; PCA selects 2000 variable genes and 20 components using the existing fixed fit plan. Mixing evaluated within source cell type and condition. Every donor held out for classifier fitting, while integration is transductive on all cells. No prospective prediction claim.',
 rareTypeRule='Evaluate every source type; rare means <1 percent of the full source cohort, defined before correction. Per-type recall retained. Missing donor/type/condition strata explicitly reported, never filled or silently pooled.',
 controlRule='Condition-centered and cell-type-centered PCA are label-informed erasure controls, never candidate corrections.',
 rejectionRule='Full Baron donor-vs-disease design is disconnected and must be rejected; never subset to force eligibility.',
 status='predeclared-before-integration-and-reference-outputs')
save('protocol.json',protocol);print(json.dumps(dict(cells=x.n_obs,features=x.n_vars,nonzeros=counts.nnz,preparedSHA256=sha(a.out/'kang.h5ad'))))
