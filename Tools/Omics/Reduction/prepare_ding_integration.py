#!/usr/bin/env python3
"""Join official public annotations to original GEO counts without inferring labels."""
import argparse,csv,gzip,hashlib,json,re
from collections import Counter
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
p=argparse.ArgumentParser(description=__doc__)
for name in ['source','annotations','base-fit','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(name,value):(a.out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
source_hashes={'cells.tsv.gz':'1011261d57c4a8c694e6c6ddebfd7107603d3789f6d5bc7ce30d91f8ecea3900','genes.tsv.gz':'ddb188e8fd1331d0098067273a882ed98bac7679a763701d838e6ff08898c00f','counts.mtx.gz':'e057891069195275e89f00c6b852380e9515367122fec06a0e563942901f09cb'}
for name,digest in source_hashes.items():assert sha(a.source/name)==digest
cells=gzip.open(a.source/'cells.tsv.gz','rt').read().splitlines();genes=gzip.open(a.source/'genes.tsv.gz','rt').read().splitlines()
assert len(cells)==len(set(cells))==44615 and len(genes)==len(set(genes))==33694
fields={}
for field in ['CellType','Experiment','Method']:
 with (a.annotations/(field+'.tsv')).open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 assert len(rows)==31021 and len({r['NAME'] for r in rows})==len(rows)
 fields[field]={r['NAME']:r[field] for r in rows}
assert fields['CellType'].keys()==fields['Experiment'].keys()==fields['Method'].keys()
# Explicit source naming aliases. DNA barcode text is retained exactly; no fuzzy,
# positional or expression-based matching. Plate/array tags stay in each key.
methods={'10x-Chromium-v2-A':'10x Chromium (v2) A','10x-Chromium-v2-B':'10x Chromium (v2) B','10x-Chromium-v3':'10x Chromium (v3)','10x-Chromium-v2':'10x Chromium (v2)','CEL-Seq2':'CEL-Seq2','Drop-seq':'Drop-seq','Seq-Well':'Seq-Well','inDrops':'inDrops','Smart-seq2':'Smart-seq2'}
prefixes={('pbmc1','10x Chromium (v2) A'):'pbmc1_10x_v2_A_',('pbmc1','10x Chromium (v2) B'):'pbmc1_10x_v2_B_',('pbmc1','10x Chromium (v3)'):'pbmc1_10x_v3_',('pbmc2','10x Chromium (v2)'):'pbmc2_10X_V2_'}
for experiment in ['pbmc1','pbmc2']:
 for method,tag in [('CEL-Seq2','Celseq2'),('Drop-seq','Drop'),('Seq-Well','Seqwell'),('inDrops','inDrops')]:prefixes[experiment,method]=experiment+'_'+tag+'_'
def source_key(name):
 experiment,method,barcode=name.split('.',2)
 return experiment.lower(),methods[method],barcode
def portal_key(name):
 experiment=fields['Experiment'][name];method=fields['Method'][name]
 if method=='Smart-seq2':return None # Read-count input is deliberately separate.
 prefix=prefixes[experiment,method];assert name.startswith(prefix)
 barcode=name[len(prefix):]
 if method=='CEL-Seq2':
  match=re.fullmatch(r'(\d+)_([ACGT]+)',barcode);assert match;barcode=match[1]+'-'+match[2]
 elif method=='Seq-Well' and experiment=='pbmc1':
  match=re.fullmatch(r'([ACGT]+)_(R\d+)',barcode);assert match;barcode=match[2]+'-'+match[1]
 elif method=='Seq-Well' and experiment=='pbmc2':
  match=re.fullmatch(r'S3_N(\d+)([ACGT]+)',barcode);assert match;barcode='S3_'+match[1]+'-'+match[2]
 elif method=='inDrops' and experiment=='pbmc2':
  match=re.fullmatch(r'1_([ACGT]+)\.([ACGT]+)\.([ACGT]+)',barcode);assert match;barcode='-'.join(match.groups())
 return experiment,method,barcode
source_keys=[source_key(c) for c in cells];assert len(set(source_keys))==len(cells)
source_lookup={k:i for i,k in enumerate(source_keys)};joined={};unmatched=[];read_annotations=[]
for name,label in fields['CellType'].items():
 key=portal_key(name)
 if key is None:read_annotations.append(name);continue
 if key not in source_lookup:unmatched.append(dict(portalCell=name,experiment=key[0],method=key[1],barcode=key[2],cellType=label));continue
 i=source_lookup[key];assert i not in joined;joined[i]=name
keep=np.array([key[1]!='Smart-seq2' for key in source_keys]);assert int(keep.sum())==44031
records=[]
for i in np.flatnonzero(keep):
 experiment,method,barcode=source_keys[i];name=joined.get(i);label=fields['CellType'][name] if name is not None else ''
 records.append(dict(sourceRow=int(i),sourceCell=cells[i],native_sample='.'.join(cells[i].split('.')[:2]),experiment=experiment,method=method,portalCell=name or '',sourceCellType=label,annotationStatus='unavailable' if name is None else ('source-unassigned' if label=='Unassigned' else 'source-assigned')))
frame=pd.DataFrame(records).set_index('sourceCell');frame.to_csv(a.out/'cell-alignment.tsv',sep='\t')
save('unmatched-portal-cells.json',unmatched);save('read-count-portal-cells.json',read_annotations)
protocol=dict(schemaVersion=1,status='declared-before-Ding-PCA-native-or-reference-results',sourceStudy='GSE132044 / SCP424',cells=44031,features=33694,seeds=[7,19,41],
 integration=dict(covariate='batch',clusters=100,diversity=2,temperature=.1,maximumIterations=10,relativeTolerance=.01,maximumWork=2000000000),modes=dict(fixed=dict(ridge=1),adaptive=dict(ridge=.2,ridgeScaling='expectedClusterBatchMass')),
 reference=dict(package='harmonypy',version='2.0.0',fixedLambda=1,adaptiveLambda=None,alpha=.2),evaluationNeighbors=30,
 gateMargins=dict(maximumCellTypeBalancedAccuracyLoss=.02,maximumPerTypeRecallLoss=.05,maximumProgramSpearmanLoss=.05,minimumNegativeControlAccuracyLoss=.10),
 programs={'NK-associated':['NKG7','GNLY','KLRD1'],'T-receptor':['CD3D','CD3E','TRAC']},programProvenance='Analyst-declared diagnostic gene sets, equal-weight mean log1p counts per 10000. Not new inferred cell labels, exclusive lineage markers, treatment responses or causal validation.',
 sourcePolicy='Every original UMI-method cell and gene enters PCA/correction. All 584 Smart-seq2 read-count cells remain an explicit incompatible-unit source. Source labels are joined by experiment, method and exact barcode including plate/array tag. No expression-based, fuzzy or row-position matches.',
 evaluationPolicy='Cell-type and classifier metrics use every available source-assigned label; unavailable and source-Unassigned cells remain in correction and full-neighbor program evaluation. Report coverage by experiment and method. Hold out each source experiment for classifier training, not an assumed donor. Mixing is evaluated within type/experiment strata. No treatment classification claim.',
 metadataPolicy='Native batchID retains the source Method. Native condition retains the source Experiment as a protected design stratum, not a treatment/health claim. biologicalReplicateID is the experiment (the paper reports one biological sample per experiment). Donor identity remains unreported. Methods shared across experiments connect this design.',
 controls='Source-cell-type centering is a label-informed diagnostic erasure control only. No thresholds or parameters selected after integration outputs.',
 qualification='Independent study evaluation of previously fixed candidate settings. Transductive correction, source-derived evaluation labels and partial annotation coverage do not establish prospective or complete biological validation. Retain prior Kang failures.')
save('protocol.json',protocol)
counts=mmread(a.source/'counts.mtx.gz').T.tocsr();assert counts.shape==(44615,33694)
assert np.issubdtype(counts.dtype,np.integer) and (counts.data>=0).all()
counts.sum_duplicates();counts.sort_indices();counts.eliminate_zeros();umi=counts[keep].tocsr()
assert (np.asarray(umi.sum(axis=1)).ravel()>0).all()
symbols=[g.partition('_')[2] for g in genes];assert all(symbols)
for program in protocol['programs'].values():
 for symbol in program:assert symbols.count(symbol)==1,symbol
var=pd.DataFrame({'symbol':symbols},index=pd.Index(genes,name='feature_id'))
ad.AnnData(X=umi,obs=frame,var=var).write_h5ad(a.out/'ding.h5ad',compression='gzip')
loaded=ad.read_h5ad(a.out/'ding.h5ad');assert (loaded.X!=umi).nnz==0 and loaded.obs_names.tolist()==frame.index.tolist() and loaded.var_names.tolist()==genes
fit=json.loads(a.base_fit.read_text());samples=[]
for sample,group in frame.groupby('native_sample',sort=True):
 assert group.experiment.nunique()==group.method.nunique()==1
 experiment=group.experiment.iloc[0];samples.append(dict(id=sample,biologicalReplicateID=experiment,condition=experiment,batchID=group.method.iloc[0],organism='NCBITaxon:9606'))
fit['featureNamespace']='Ding2019/GRCh38/original-Ensembl-symbol-IDs'
fit['mapping']=dict(schemaVersion=1,id='ding2019-all-umi-cells',evidence='measured',countUnit='umiCount',matrixPath='X',sampleColumn='native_sample',featureNameColumn='symbol',mitochondrialFeatureIDs=[],samples=samples,sourceDescription='Ding et al. GSE132044: all 44031 original UMI-method cells and all 33694 genes. Smart-seq2 read counts excluded by unit, never relabeled UMI. Source Experiment is the protected native condition/design stratum; Method is batch covariate. Donor identity unreported. CellType annotations are evaluation-only and not mapped into native groups.')
save('fit.json',fit)
coverage=frame.groupby(['experiment','method','annotationStatus'],observed=True).size().reset_index(name='cells').to_dict('records')
save('checks.json',dict(status='passed',sourceFiles=source_hashes,annotationFiles={field:sha(a.annotations/(field+'.tsv')) for field in fields},preparedSHA256=sha(a.out/'ding.h5ad'),originalCells=44615,umiCells=44031,readCountCells=584,features=33694,originalNonzeros=int(counts.nnz),umiNonzeros=int(umi.nnz),allRetainedCountValuesExact=True,allRetainedCellFeatureIdentitiesExact=True,sourceLabelMatches=len(joined),unmatchedPortalUMICells=len(unmatched),unmappedReadCountAnnotations=len(read_annotations),coverage=coverage,sourceCellTypeCounts=frame.sourceCellType.value_counts().to_dict(),qualification=protocol['qualification']))
print((a.out/'checks.json').read_text())
