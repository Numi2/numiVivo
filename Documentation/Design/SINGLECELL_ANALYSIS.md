# Native single-cell cohort analysis

## Implemented route

The count workflow now accepts bounded plain/gzip MEX source libraries and keeps
the original source bytes in the existing immutable artifact store. A separate
analysis plan applies explicit cell selection, retains source-row lineage,
computes feature summaries, creates the existing raw pseudobulk representation
and evaluates named condition contrasts. The public CLI can reconstruct and
export the full analysis, selected raw MEX libraries and inspectable TSV tables.
No Python, R, CUDA, shell decompressor or external statistics runtime executes
production calculations. Native gzip uses the platform zlib library.

`VivoSingleCellProcessing` owns cell decisions and sparse feature summaries.
`VivoPseudobulkDifferentialExpression` builds statistical units, size factors and
designs. `VivoOmicsQR`/`VivoOmicsLinearStatistics` own the bounded FP64 numerical
model. `VivoSingleCellCohortAnalysis` composes these authorities. Dedicated CLI
receipts and `vivo.platform.singlecell-analyze` both reuse this composition and
the pre-existing store/scheduler; neither creates a separate data platform.

## Cell and feature semantics

Filters have explicit minimum/maximum counts and detected-feature bounds,
optional mitochondrial-fraction bounds, a declared missing-mitochondrial policy,
an optional sample subset and explicitly excluded sample/barcode identities.
Every source cell receives an acceptance decision, rejection reasons and its
original quality metrics. A missing annotation is not zero mitochondrial content.
The default filter only requires at least one count and one detected feature;
it is not an assay-specific quality recommendation or a doublet detector.

All original features remain in the selected raw matrix. Sparse feature summaries
report exact total counts and detected cells, plus means/variances of the separate
library-size/log1p view. Implicit zeros are included using a Welford accumulator
combined with the zero population. No dense genes-by-cells matrix is formed.
Empty selections have missing means/variances. Pseudobulk source indices refer
to selected rows and are mapped back using the retained selected-to-source map.

## Statistical model

Each named contrast selects one supplied cell-group label and two conditions.
The default minimum is three biological replicates per condition and ten cells
per pseudobulk. These are admission requirements, not power calculations. Group
labels, donor identities and independence are supplied metadata, not discoveries.

Independent designs reject repeated biological-replicate or known donor IDs.
Paired designs require exactly one usable pseudobulk per donor and condition and
include donor fixed effects. Incomplete pairs are errors. Optional donor subsets
are explicit. Categorical batch columns are included by default; a pseudobulk
spanning multiple batches cannot receive a fabricated batch label. A rank-deficient
or overparameterized design is rejected instead of dropping a confounded column.
Turning batch adjustment off is an explicit modelling decision, retained in the
request; it does not establish the absence of batch effects.

Library factors are either declared total-library factors or median ratios over
genes positive in every selected pseudobulk. Median-ratio factors use arithmetic
medians of ratios, including even reference-set sizes, then geometric-mean
rescaling. The selected reference features, raw library counts and final factors
are recorded. Too few reference genes is an error, not automatic normalization
fallback. The expression response is `log2(count / factor + priorCount)`.
Counts and library totals above 2^53 are rejected for inference; raw import and
export still preserve UInt64 values exactly.

The linear model uses twice-reorthogonalized QR, not normal equations. The
condition contrast, full design, coefficient names and statistical observations
are retained. The optional empirical-Bayes method fits a common scaled
inverse-chi-square variance prior by log-variance moments. It uses at least 20
positive residual variances; near-zero variances (<=1e-24) are excluded from prior
estimation. Prior degrees of freedom are bounded at 1,000,000, with the boundary
reported. Posterior variance combines residual and prior sums of squares.
Testing degrees of freedom are capped by the pooled residual degrees of freedom.

Ordinary or moderated Student t inference reports the fitted log2 expression
effect, standard error, t statistic, degrees of freedom, interval and p-value.
Filtered features have no invented test result. An unmoderated zero-variance fit
has an effect but no fabricated p-value. Benjamini-Hochberg correction is applied
only across tested features within each named contrast. It does not provide a
joint across-contrast guarantee or repair selection after inspecting results.

This is an untrended, unweighted log-expression model. It does not implement voom
mean-variance weights, robust/trended variance moderation, negative-binomial
likelihoods, random-effects models or validated DESeq2/edgeR/limma equivalence.
Its inferential assumptions require experimental evaluation; synthetic signal
recovery and correct arithmetic do not establish false-discovery calibration.

## Interchange, lineage and replay

Gzip decoding validates member checksums/trailers, supports up to 1,024
concatenated members, rejects trailing non-gzip data and enforces an aggregate
expanded-byte allowance independently of the aggregate original-byte allowance.
Both defaults are 64 MiB across the entire campaign, including the manifest.
Cell, feature, nonzero and line limits remain active. This is bounded in-memory
processing, not out-of-core or atlas-scale execution.

MEX export preserves raw integer counts, ordered feature identities, sample and
barcode identities, supplied groups and mitochondrial annotations. Files are
organized by sample, and a sidecar preserves exported-row-to-source-row mapping
when original rows were interleaved. Analysis exports also retain the selected-
to-original map. New-directory exports never overwrite an existing destination.
TSV exports carry a receipt index and use empty fields for missing values.

Analysis receipts bind the existing source count receipt, exact analysis-plan
bytes, the result and the executable/OS implementation identity. Verification
checks hashes and recursively reconstructs counts and analysis. A correctly
hashed but altered statistical report is not accepted. Source byte changes
create a different input identity; existing snapshots remain intact. The common
workflow adapter validates outputs through deterministic reconstruction and
retains existing cache and dependency semantics.

Current limits include 32 contrasts, 512 observations and 128 design columns per
contrast, 250 million design-response products per contrast, and explicit source,
plan, artifact and export byte caps. These are allocation/work admission limits,
not measured performance or resident-memory guarantees.

## Execution scope and remaining work

The standalone new FP64 statistics and system-zlib decoder were compiled and
smoke-executed with Swift 6.2.1 on x86_64 Linux during development. The remaining
new Swift integration files were syntax-parsed. The complete Apple product and
new end-to-end CLI route were not executed in that development environment.
`Single-cell count workflows` contains the Apple cohort suites and real-CLI
checker, but adding that gate is not evidence that its latest run has passed.
The [scoped execution record](../Audit/2026-09-08_SINGLECELL_COHORT.json) records
the standalone command, source hashes and observed output.

Not implemented here: HDF5/AnnData, multimodal feature selection, gene-ID
harmonization, out-of-core arrays, doublet/ambient correction, clustering, inferred
cell types, GPU acceleration or externally calibrated biological predictions.
The molecular/prepared-reaction release work is unaffected by this increment.

## Method and format references

- 10x Genomics MEX format: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-outputs-mex-matrices
- zlib gzip/concatenated-member semantics: https://www.zlib.net/manual.html
- Squair et al., biological replication in single-cell differential expression: https://www.nature.com/articles/s41467-021-25960-2
- Scaled-F variance-prior methodology reference (not runtime dependency): https://github.com/bioc/limma/blob/devel/R/fitFDist.R
