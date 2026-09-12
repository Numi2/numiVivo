# Paired reconstruction diagnostics: neither gate passes

The ATAC reconstruction diagnostic also fails its combined development
no-worsening rule. Together with the retained [RNA result](../RNAPreservation/README.md),
this withholds biological-preservation promotion for the current weighted graph.
Numerical and source-replay qualification remain valid.

| Weighted versus modality-only baseline | Primary neighbor-mean SSE change | Secondary fuzzy-mean SSE change | Combined gate |
| --- | ---: | ---: | --- |
| RNA | +0.8527% | -0.1274% | FAIL |
| ATAC | +0.04437% | -0.30655% | FAIL |

Positive changes mean worse reconstruction. Both primary comparisons worsen
across all three detection-frequency groups. Both secondary comparisons improve,
but the frozen rules require no worsening for both uses; favorable secondary
results cannot replace the primary results.

## Complete ATAC evaluation

All 2,711 nuclei and all 98,319 peaks are evaluated, including the one zero peak.
The source contains 19,292,713 nonzero measurements and 48,245,242 cut sites.
Native product cell/feature totals are independently checked against the source
before count-based method-1 TF-IDF is computed. No binarization, peak selection,
cell selection or conversion from cut sites to fragments is performed.

Targets are the complete TF-IDF values. Predictions use either the mean of
14 non-self selected neighbors or row-normalized native fuzzy connectivity.
The four fixed graphs are unchanged from earlier qualification. Sparse 128-feature
blocks avoid constructing a complete dense cells-by-peaks matrix.

| Graph | Neighbor-mean SSE | Fuzzy-mean SSE |
| --- | ---: | ---: |
| RNA-only | 58,653,472.64 | 59,863,538.96 |
| ATAC-only | 57,898,358.99 | 59,065,163.53 |
| Equal-weight joint | 57,991,406.13 | 59,134,238.52 |
| Cell-specific weighted | 57,924,047.98 | 58,884,099.72 |
| Leave-one-cell-out global feature mean | 56,887,850.59 | 56,887,850.59 |

Weighted reconstruction is worse than the simple leave-one-out mean by 1.82%
and 3.51% in aggregate, respectively. Macro variance-normalized errors are worse
by 2.09% and 3.80%. This does not demonstrate useful general prediction of peak
measurements, even though some comparisons among graphs improve.

There are 5,738 peaks detected in 1–19 cells, 53,310 in 20–99 cells and 39,270
in at least 100 cells. Macro error includes 92,580 variable peaks detected in
at least 20 cells. The zero peak is retained with undefined normalized error.

All-feature squared errors agree between two algebraic calculations within the
declared floating-point checks (maximum absolute difference 7.28e-12).
Independent scalar neighbor averaging across all cells for 17 prespecified peaks
agrees within 1.48e-12. Complete per-feature arrays, baselines, null comparisons,
versions, protocol and code are retained in the verified archive/manifest.
Inherited RNA counter labels were corrected to feature labels without changing
numerical values or gates; the original report and code remain retained.

## Interpretation and next hypothesis

All peaks were available to LSI. This is an in-sample reconstruction diagnostic,
not independent validation of ATAC biology. One donor cannot establish donor
integration, unseen-context prediction, regulatory causality or clinical effects.
Measurement noise and representation choices remain possible explanations; these
results do not identify a biological mechanism.

A concrete next preprocessing hypothesis is RNA component scaling. The current
paired workflow standardizes RNA PC columns to unit variance.
[Seurat 5.3.0 RunPCA](https://github.com/satijalab/seurat/blob/v5.3.0/R/dimensional_reduction.R)
defaults to retaining singular-value weighting in cell scores. Test an explicit
variance-preserving RNA profile while retaining legacy source replay and both
modalities' fixed/matched baselines. Other preprocessing differences remain;
this observation is not a claim of full Seurat parity or an improved outcome.

The current results must remain visible in that development comparison. Any
candidate improvement on these already inspected data still requires separately
specified experimental evidence before broad biological promotion.

The [variance-retaining candidate](../RNAVarianceProfile/README.md) has now been executed. Both modality gates still fail; its full results and numerical-verification limits are retained.
