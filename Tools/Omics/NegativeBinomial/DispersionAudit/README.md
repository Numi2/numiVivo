# Dispersion-stage diagnosis of the measured null result

Follow-up: the [profile audit](../ProfileAudit/README.md) finds missed better
objectives in the default reference and passes the selected native numerical
checks. The differences below should not be read as proof of a native numerical
defect. Reference optimization sensitivity and calibration remain distinct.

This post-score audit follows the [frozen untreated-cell benchmark](../NullBenchmark/README.md).
It narrows the next inference work to gene-wise dispersion estimation, its
reference-gene admission and outlier handling. It does not change or qualify
the native model, its original calls, or general false-discovery control.

All twenty pinned DESeq2 refits reproduce every archived final dispersion,
effect, standard error, p-value and adjusted p-value exactly after parsing:
maximum absolute difference **0**. Outlier and convergence flags also match.
All eighty declared trend fits complete without warnings or failed cases.
Native inputs come from the full model archives, with both compressed and
logical hashes checked on the owning Mac mini before extracting stage tables.

| Diagnostic across ten splits | Kang | Hagai |
|---|---:|---:|
| Native prior log variance | 0.891–1.037 | 0.25 throughout |
| DESeq2 prior log variance | 0.690–0.960 | 0.25 throughout |
| Native residual log variance | 1.222–1.367 | 1.368–1.556 |
| DESeq2 residual log variance | 1.021–1.290 | 0.360–0.398 |
| Median DESeq2/native trend dispersion per split | 2.017–2.549 | 2.536–2.857 |
| Median DESeq2/native final dispersion per split | 1.582–2.044 | 2.365–2.640 |

The dispersion ratios use genes with both stage estimates and retain each
method's coverage separately. Many low-boundary gene-wise estimates agree;
medians alone obscure the substantial disagreement in interior membership.
For example, Hagai seed 1 has 4,685 native trend reference genes, while DESeq2
has 1,258 estimates above its 1e-6 admission floor. Among jointly reported genes,
3,506 are below that floor in DESeq2 and above it in native.

## What the decomposition establishes

The [stage protocol](PROTOCOL.md) and separately declared
[trend decomposition](TREND_PROTOCOL.md) preserve the original data and inference.
Pinned DESeq2's trend routine is applied to four fixed input sets in every split.
Its original-input case reproduces its stored curve. On native gene-wise
estimates and native reference genes, it nearly reproduces the native Hagai
curve: maximum relative difference over observed means is 2.31e-5 across all
ten splits. The corresponding Kang maximum is 2.70%; that smaller algorithmic
discrepancy remains recorded, not rounded into an exact agreement claim.

Replacing only the gene-wise estimates with DESeq2 estimates on the same native
reference set raises the median Hagai curve to 2.428–2.734 times the native curve
and to 0.913–0.987 of the original DESeq2 curve. Conversely, merely applying
DESeq2's 1e-6 admission floor to native estimates leaves the median curve within
0.043% of its original value. Thus the threshold alone does not explain the
large curve difference. The gene-wise estimation procedure is a concrete next
target for an independently checked likelihood/profile comparison.

Both methods use a 0.25 prior variance on every Hagai split, so changing prior
strength cannot explain their current disagreement. Their residual-spread
estimates used for outlier admission differ substantially. DESeq2 flags five
of the twenty observed native Hagai calls as dispersion outliers; native flags
none. The remaining fifteen calls show that outlier handling alone is not a
complete explanation either. Whole-cohort and per-call stage tables are retained.

Native currently refits coefficients at each Cox–Reid profile point, whereas
DESeq2 estimates dispersions using previously fitted means; the procedures are
not identical estimators. These measured differences do not prove either one
correct, justify copying a threshold to reduce calls, or replace an independent
calibration/power experiment after a principled method change.

## Reproduce and inspect

Use the same qualified Python environment and R library as the original null
benchmark. The runner reads its retained input tables and checks the original
DESeq2 outputs. It needs SSH access to the hash-bound native model archives.
The default Rscript path and native archive host are explicit in `run.py`.

```sh
python run.py --root /runs/null-benchmark --out /new/stage-audit \
  --r-library /qualified/R-library
R_LIBS_USER=/qualified/R-library Rscript trends.R \
  /new/stage-audit /new/stage-audit/trend-decomposition
python summarize.py --root /new/stage-audit
```

The [evidence manifest](evidence/2026-09-10/manifest.json) binds all twenty stage
tables, native extracts, reference session information, per-fit logs, eighty
curve outputs, per-call diagnostics and summaries. The original null benchmark
retains count matrices, cell assignments and full model archive hashes. No new
Swift build is claimed: this audit reads the already qualified native binary's
results at 6e6f8b2 and uses Python/R only as reference and diagnostic tooling.
