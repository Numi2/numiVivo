#!/usr/bin/env python3
"""Prepare the first measured RNA assay; requires anndata, numpy and pandas.

Input is the deposited Kang benchmark H5AD, never simulated expression.
Preparation preserves the full measured gene axis and original source hash.
"""
import argparse
from pathlib import Path
import shutil

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from wetlab import FORMAT, write, sha, native, read, require, catalog

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('source', type=Path)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
source = ad.read_h5ad(a.source)
require(source.shape == (24673, 15706), 'Expected the documented Kang 24,673-cell / 15,706-gene benchmark')
require(sha(a.source) == 'e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830',
        'Kang input differs from the source-bound benchmark; review provenance before adding a new source')
selected = np.flatnonzero(source.obs.cell_type.astype(str).to_numpy() == 'B cells')
x = source.X[selected].tocsr()
require(np.all(np.isfinite(x.data)) and np.all(x.data >= 0) and np.all(x.data == np.floor(x.data)),
        'Expected measured integer UMI counts')
obs = source.obs.iloc[selected]
donors = obs.replicate.astype(str); conditions = obs.label.astype(str)
require(set(conditions) == {'ctrl', 'stim'}, 'Unexpected conditions')
require(len(set(donors)) == 8, 'Expected eight source-labelled donors')
sample_ids = donors + '__' + conditions
samples = [{'id': d + '__' + c, 'donorID': d, 'biologicalReplicateID': d,
            'condition': c, 'batchID': 'unreported', 'organism': 'NCBITaxon:9606'}
           for d in sorted(set(donors)) for c in ('ctrl', 'stim')]
prepared = ad.AnnData(X=csr_matrix(x, dtype=np.uint32),
    obs=pd.DataFrame({'sample': sample_ids.to_numpy(), 'cellGroup': 'B cells'}, index=obs.index),
    var=pd.DataFrame(index=source.var_names))
prepared.write_h5ad(a.output / 'b-cells.h5ad', compression='gzip')
write(a.output / 'preparation.json', {
    'originalSHA256': sha(a.source), 'preparedSHA256': sha(a.output / 'b-cells.h5ad'),
    'originalObservationIndices': selected.tolist(), 'cells': len(selected), 'features': source.n_vars,
    'countsPreserved': True, 'sourceAnnotations': ['cell_type', 'replicate', 'label'],
    'preparationScriptSHA256': sha(__file__),
    'limits': 'Deposited donor labels; primary library accession mapping unresolved; batch unreported.'})
# Retain the source for auditing the annotation-only preparation.
shutil.copy2(a.source, a.output / 'original-kang.h5ad')
mapping = {'schemaVersion': 1, 'id': 'kang-b-cells', 'evidence': 'measured',
           'sourceDescription': 'Kang 2018 GSE96583 deposited B-cell annotations; eight paired donor labels. '
                                'Study-level 6h exposure; no dose interpolation. Batch unreported.',
           'countUnit': 'umiCount', 'matrixPath': 'X', 'samples': samples,
           'sampleColumn': 'sample', 'groupColumn': 'cellGroup'}
write(a.output / 'source-plan.json', {'schemaVersion': 1, 'mapping': mapping, 'contrasts': []})
native(a.binary, ['singlecell-h5ad-pseudobulk', a.output / 'b-cells.h5ad', '--plan',
       a.output / 'source-plan.json', '--output', a.output / 'source'], a.output / 'prepare-log.json')
write(a.output / 'assay.json', {
    'format': FORMAT, 'id': 'kang-b-cell-ifnb-6h', 'title': 'Immune-cell RNA response',
    'family': 'cell-response', 'cellGroup': 'B cells', 'intervention': 'IFN-beta',
    'observationHours': 6, 'controlCondition': 'ctrl', 'treatmentCondition': 'stim',
    'featureNamespace': 'HUMAN_GENE_SYMBOL', 'minimumCellsPerArm': 20,
    'sourceBundle': 'source', 'sourceReportSHA256': sha(a.output / 'source/report.json'),
    'preparationSHA256': sha(a.output / 'preparation.json'),
    'sourceCitation': 'https://stacks.cdc.gov/view/cdc/79371/cdc_79371_DS1.pdf',
    'provenance': 'Kang 2018 GSE96583 public benchmark, source-labelled B cells from eight lupus donors. '
                  'Previously inspected development cohort. Study-level six-hour exposure; source library '
                  'accession mapping remains unresolved and batch is unreported. Original H5AD SHA256: ' + sha(a.source)})
print(catalog(a.output / 'assay.json'))
