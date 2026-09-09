#!/usr/bin/env python3
"""Real Kang 2018 B-cell paired-donor benchmark: native baseline vs Scanpy/PyDESeq2.

All B cells and all source genes are retained. Counts are never reconstructed
from normalized expression. Cells are not treated as independent DE replicates.
"""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess
import time
import warnings
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

SOURCE_URL = 'https://ndownloader.figshare.com/files/34464122'
SOURCE_SHA = 'e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830'
EXPECTED = ['ISG15','IFIT1','IFIT3','MX1','OAS1']
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
commands = []

def save(name, obj):
    path = a.out / name
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    return path

def run(label, *arguments):
    command = [str(a.binary.resolve()), *map(str, arguments)]
    start = time.perf_counter()
    result = subprocess.run(['/usr/bin/time', '-l', *command], capture_output=True, text=True)
    (a.out / (label + '.log')).write_text(result.stdout + result.stderr)
    peak = next((int(line.split()[0]) for line in result.stderr.splitlines() if 'maximum resident set size' in line), None)
    commands.append(dict(label=label, command=command, returncode=result.returncode, elapsedSeconds=time.perf_counter()-start, peakResidentBytes=peak))
    save('commands.json', commands)
    if result.returncode:
        save('failure.json', dict(status='native-command-failed', command=commands[-1], scope='Kang2018 all B cells; no synthetic fallback'))
        raise RuntimeError(label + ': ' + result.stderr)

assert hashlib.sha256(a.source.read_bytes()).hexdigest() == SOURCE_SHA, 'wrong source release'
upstream = ad.read_h5ad(a.source, backed='r')
assert upstream.shape == (24673,15706)
rows = np.flatnonzero(upstream.obs['cell_type'].astype(str).to_numpy() == 'B cells')
counts = sparse.csr_matrix(upstream[rows,:].X)
assert counts.shape == (2651,15706) and counts.nnz == 1479543
assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and np.equal(counts.data, np.floor(counts.data)).all()
counts = counts.astype(np.int64)
obs, var = upstream.obs.iloc[rows].copy(), upstream.var.copy()
upstream.file.close()
assert var.index.is_unique and obs.index.is_unique
obs['native_sample'] = obs['replicate'].astype(str) + '__' + obs['label'].astype(str)
assert sorted(obs['label'].unique()) == ['ctrl','stim'] and obs['replicate'].nunique() == 8
pairs = sorted(set(zip(obs['replicate'].astype(str), obs['label'].astype(str))))
assert len(pairs) == 16
samples = [dict(id=donor+'__'+condition, biologicalReplicateID=donor, donorID=donor, condition=condition,
                batchID='unreported', organism='NCBITaxon:9606') for donor,condition in pairs]
prepared = ad.AnnData(X=counts, obs=obs, var=var)
prepared.uns['numivivo_benchmark_source'] = dict(url=SOURCE_URL, sha256=SOURCE_SHA, row_indices=rows,
    scope='All B cells, all genes; no expression-based cell or feature subsampling', donor_key='replicate', condition_key='label')
prepared_path = a.out / 'prepared.h5ad'
prepared.write_h5ad(prepared_path, compression='gzip')
manifest = dict(schemaVersion=1, id='kang2018-bcells', evidence='measured', countUnit='umiCount', matrixPath='X',
    sourceDescription='Kang 2018 (doi:10.1038/nbt.4042), pertpy file 34464122; all 2651 B cells and all 15706 source genes. Batch metadata unreported; paired donor model does not adjust for batch.',
    sampleColumn='native_sample', groupColumn='cell_type', samples=samples, mitochondrialFeatureIDs=[])
map_path = save('mapping.json', manifest)
analysis_path = save('analysis.json', dict(schemaVersion=1, id='kang2018-bcells-paired', normalizationTarget=10000,
    filter=dict(minimumCounts=1, minimumDetectedFeatures=1), contrasts=[dict(id='IFNB-vs-control', controlCondition='ctrl',
    treatmentCondition='stim', cellGroup='B cells', design='pairedDonors', sizeFactors='medianRatio', variance='empiricalBayes',
    adjustForBatch=False, minimumCellsPerPseudobulk=10, minimumReplicatesPerCondition=3, minimumFeatureCounts=10,
    minimumExpressingPseudobulks=3, minimumReferenceFeatures=10, priorCount=0.5)]))
save('design-and-acceptance.json', dict(schemaVersion=1, sourceSHA256=SOURCE_SHA, publication='https://doi.org/10.1038/nbt.4042',
    cellType='B cells', donors=8, design='~ donor + condition', batchAdjustment=False, batchMetadata='not provided in input',
    cellsPerSample={d+'__'+c:int(((obs.replicate==d)&(obs.label==c)).sum()) for d,c in pairs},
    expectedPositiveGenes=EXPECTED, expectedDirectionRequirement='at least 4 of 5 genes have positive log2 effect in both models',
    comparisonThreshold='No parity/competitiveness threshold; different count models are compared descriptively',
    referencePolicy=dict(refitCooks=False, cooksFilter=False, independentFiltering=False,
        reason='Controlled same-count, same-filter comparison; robust-reference sensitivity is not qualified by this run')))
imported, store = a.out / 'imported', a.out / 'store'
run('native-import', 'singlecell-h5ad-import', prepared_path, '--plan', map_path, '--output', imported)
run('native-counts', 'singlecell-run', imported/'manifest.json', '--store', store, '--output', a.out/'count-receipt.json')
run('native-analysis', 'singlecell-analyze', a.out/'count-receipt.json', '--plan', analysis_path, '--store', store, '--output', a.out/'analysis-receipt.json')
run('native-replay-export', 'singlecell-analysis-export', a.out/'analysis-receipt.json', '--store', store, '--output', a.out/'native-report.json')
native = json.loads((a.out/'native-report.json').read_text())
processed = native['processed']
# MEX publication groups rows by sample. Compare by biological identities,
# never assume that the input file and its MEX projection share row ordering.
identity = {str(barcode):i for i,barcode in enumerate(obs.index)}
native_order = np.array([identity[cell['barcode']] for cell in processed['dataset']['cells']])
assert len(native_order) == len(rows) and len(set(native_order)) == len(rows)
nm = processed['dataset']['matrix']
recovered = sparse.csr_matrix((np.array(nm['counts'],dtype=np.int64), nm['featureIndices'],nm['rowOffsets']),shape=counts.shape)
assert (recovered != counts[native_order]).nnz == 0
ref = ad.AnnData(X=counts.astype(np.float64).copy(), obs=obs.copy(), var=var.copy())
sc.pp.calculate_qc_metrics(ref, percent_top=None, log1p=False, inplace=True)
qc = [decision['quality'] for decision in processed['decisions']]
np.testing.assert_array_equal([q['totalCounts'] for q in qc], ref.obs.total_counts.to_numpy()[native_order])
np.testing.assert_array_equal([q['detectedFeatures'] for q in qc], ref.obs.n_genes_by_counts.to_numpy()[native_order])
sc.pp.normalize_total(ref, target_sum=10000)
sc.pp.log1p(ref)
normal = processed['normalized']
native_normal = sparse.csr_matrix((normal['values'], normal['featureIndices'],normal['rowOffsets']), shape=counts.shape)
ref_normal = ref.X[native_order].tocsr()
assert np.array_equal(native_normal.indptr,ref_normal.indptr) and np.array_equal(native_normal.indices,ref_normal.indices)
normal_error = float(np.max(np.abs(native_normal.data-ref_normal.data)))
assert normal_error < 1e-10
bulk_rows = [np.asarray(counts[((obs.replicate==donor)&(obs.label==condition)).to_numpy()].sum(axis=0)).ravel() for donor,condition in pairs]
# Only the small donor-by-gene aggregate is dense; cells-by-genes stays sparse.
bulk = np.vstack(bulk_rows).astype(np.int64)
nb = processed['pseudobulk']
assert [(g['donorID'],g['condition']) for g in nb['groups']] == pairs
bm = nb['matrix']
native_bulk = sparse.csr_matrix((np.array(bm['counts'],dtype=np.int64),bm['featureIndices'],bm['rowOffsets']),shape=bulk.shape)
np.testing.assert_array_equal(native_bulk.toarray(),bulk)
eligible = (bulk.sum(axis=0)>=10) & ((bulk>0).sum(axis=0)>=3)
metadata = pd.DataFrame(dict(donor=[d for d,c in pairs],condition=[c for d,c in pairs]),index=[d+'__'+c for d,c in pairs]).astype('category')
reference_counts = pd.DataFrame(bulk[:,eligible],index=metadata.index,columns=var.index[eligible])
reference_counts.to_csv(a.out/'reference-pseudobulk.tsv',sep='\t')
metadata.to_csv(a.out/'reference-design.tsv',sep='\t')
start = time.perf_counter()
dds = DeseqDataSet(counts=reference_counts, metadata=metadata, design='~ donor + condition', refit_cooks=False, n_cpus=2, quiet=True)
with warnings.catch_warnings(record=True) as reference_warnings:
    warnings.simplefilter('always')
    dds.deseq2()
    ds = DeseqStats(dds, contrast=['condition','stim','ctrl'], cooks_filter=False, independent_filter=False, n_cpus=2, quiet=True)
    ds.summary()
reference_diagnostics = [dict(category=w.category.__name__, message=str(w.message)) for w in reference_warnings]
save('reference-diagnostics.json', reference_diagnostics)
reference = ds.results_df
reference.to_csv(a.out/'pydeseq2.tsv',sep='\t')
reference_seconds = time.perf_counter()-start
native_de = pd.DataFrame(native['contrasts'][0]['features']).set_index('featureID')
native_de.to_csv(a.out/'native-de.tsv',sep='\t')
comparison = native_de.join(reference,how='inner',rsuffix='_reference')
finite = comparison[(comparison.status=='tested') & np.isfinite(comparison.log2FoldChange) & np.isfinite(comparison.log2FoldChange_reference)]
native_top = set(finite.nsmallest(50,'adjustedPValue').index)
reference_top = set(finite.nsmallest(50,'padj').index)
biology = {g:dict(nativeLog2Effect=float(native_de.loc[g,'log2FoldChange']), referenceLog2Effect=float(reference.loc[g,'log2FoldChange']),
    nativeBH=float(native_de.loc[g,'adjustedPValue']), referenceBH=float(reference.loc[g,'padj'])) for g in EXPECTED}
positive = sum(v['nativeLog2Effect']>0 and v['referenceLog2Effect']>0 for v in biology.values())
versions = {x:version(x) for x in ['scanpy','anndata','pydeseq2','numpy','scipy','pandas','h5py']}
report = dict(schemaVersion=1,status='completed-real-data-comparison',sourceSHA256=SOURCE_SHA,
    binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),versions=versions,platform=platform.platform(),
    cells=len(rows),genes=counts.shape[1],nonzeros=counts.nnz,donors=8,pseudobulks=16,
    integrity=dict(exactImportedCounts=True,exactPseudobulkCounts=True,exactScanpyQC=True,maxLogNormalizationError=normal_error),
    modelComparison=dict(native='moderated log-linear baseline, not native NB',reference='PyDESeq2 negative binomial',
        eligibleGenes=int(eligible.sum()),nativeTested=int((native_de.status=='tested').sum()),referenceFinite=int(np.isfinite(reference.pvalue).sum()),
        commonFiniteGenes=len(finite),effectSpearman=float(stats.spearmanr(finite.log2FoldChange,finite.log2FoldChange_reference).statistic),
        effectSignAgreement=float(np.mean(np.sign(finite.log2FoldChange)==np.sign(finite.log2FoldChange_reference))),
        top50BHOverlap=len(native_top&reference_top),referenceFitSeconds=reference_seconds, referenceDiagnostics=reference_diagnostics),
    expectedBiology=dict(status='passed' if positive>=4 else 'failed',positiveInBoth=positive,required=4,genes=biology),
    qualificationBoundary='One predefined B-cell IFNB contrast in one real paired-donor study; not multi-dataset competitiveness, causal proof or native NB qualification',
    nativeCommands=commands)
save('report.json',report)
print(json.dumps({k:v for k,v in report.items() if k!='nativeCommands'},indent=2))
if positive<4: raise SystemExit(1)
