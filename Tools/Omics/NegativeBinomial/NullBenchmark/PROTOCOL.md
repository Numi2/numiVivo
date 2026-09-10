# Untreated-cell split null benchmark

Declared 2026-09-10 before preparing any split or fitting its outcomes.
Native owner starts at f40b37c870e6e1702e3feff8f7b400188cce82ca.

Use the complete previously qualified Kang B-cell and Hagai mouse H5AD sources
and their original donor/condition mappings. Retain every source feature and
original observation. Only original untreated cells enter the sham contrasts:
Kang `ctrl`, Hagai `unstimulated`. Do not pool species, substitute cells for
biological replicates, infer extra donors, or combine real treatment conditions.

For each study, use exactly seeds 1 through 10. Within each original untreated
sample, sort cells by SHA256 of the UTF-8 string
`numivivo-null-v1|study|seed|sampleID|barcode`, breaking any hash tie by original
row index. Assign alternating ranks to `shamA` and `shamB`. Both arms retain
original donor, biological-replicate, cell-type and batch provenance. Split arms
are paired subsamples of the same original untreated library, not independent
new experimental samples. Each source cell belongs to exactly one arm per seed.
No counts, gene identities, annotations or expression effects enter assignment.

Use native source-preserving H5AD annotation and source-bound cell selection.
Retain complete source counts and explicit row assignments. Independently verify
all arm counts and source totals with sparse SciPy aggregation. No dense
cells-by-genes representation. Only donor/arm-by-gene aggregates may be dense.

Run the native NB contrast with paired donor fixed effects, minimum three donors,
minimum ten cells per arm, minimum total gene count ten and expression in at
least three pseudobulks. Use native median-ratio factors, gammaParametric
dispersion trend, minimum 20 trend genes, minimum prior variance 0.25 and outlier
threshold two standard deviations. Run both original full-support policy and
explicit active-donor profile policy; neither may silently replace the other.
No influence filtering, count replacement, expression-guided split rejection,
trend fallback, optimization tuning or post-result seed selection.

Compare R edgeR robust QL, limma-voom robust empirical Bayes and DESeq2 Wald using
the same paired numeric design and eligible genes. Reuse the pinned direct
Bioconductor comparison implementation and both normalization modes (native
size factors and package normalization). Preserve package warnings/failures,
DESeq2's ordinary filtering output and the explicit unfiltered comparison.
Do not change package versions or fit methods after viewing results.

For every study, seed and method report attempted/eligible/tested/withheld genes,
fit failures and boundaries, raw p-value fractions below 0.01 and 0.05, BH
rejections at 0.01/0.05/0.10, and whether any BH<0.05 rejection occurs. All tested
sham associations are false discoveries under the randomized assignment null.
Retain method-specific families and a separate jointly tested family; do not
promote intersection-only agreement. Report per-seed values and descriptive
aggregate counts without treating genes or overlapping resamples as independent.

Expected population effects are zero under the assignment randomization; the
finite realized arm difference is sampling variation. These nulls condition on
existing untreated cell pools, their sampling/capture noise and original study
selection. They do not recreate new independently sampled donors, prove biological
nulls for original treatments, establish alternative-model power, or prove general
selection-adjusted FDR control. Ten overlapping splits per study provide a stress
test, not a precise experiment-level FDR estimate. Any software fit/reference
failure remains a failure; no reduced cohort is substituted to obtain a pass.

Primary methodological context: [Squair et al.](https://www.nature.com/articles/s41467-021-25960-2)
and the [DESeq2 workflow](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).
This exact hash-split experiment is our declared diagnostic, not a reproduction
of either publication's complete evaluation protocol.
