# Explicit active-donor NB profiles

An opt-in `negativeBinomialOptions.zeroTotalDonorPolicy: "activeDonorProfile"`
handles a gene's complete zero-total donor pairs in a paired contrast. The
existing no-policy behavior is unchanged. This is an experimental inferential
method with real-data numerical checks; selection-adjusted FDR and broader
biological calibration remain open.

## Method and interpretation

A pair is removed from the gene's inferential design only when both original raw
counts are exactly zero. No individual observation, nonzero low count or partial
pair is removed. The retained pairs must satisfy the original requested minimum
replication (at least three donors). The shared design builder rebases donor and
optional batch dummy columns. Full design rank and full positive-count support
rank are still required. Independent-replicate designs reject this option.

The zero pair's unadjusted NB log likelihood and treatment Schur-complement
information vanish as its nuisance rate approaches zero. This motivates an
explicit reduced active-donor profile. It does **not** establish equality between
the singular full-design Cox–Reid determinant limit and the reduced adjusted
profile, nor establish unbiased inference after response-dependent selection.
The reported method name distinguishes this policy from the original method.

Gene-wise dispersion and the final log-prior MAP dispersion are estimated with
retained counts, original size-factor offsets and the rebuilt design. Expected
Fisher covariance, Wald inference and influence diagnostics use those same rows.
The dispersion prior is estimated from the unchanged original full-support
reference cohort, with its original residual degrees of freedom. Reduced profiles
borrow that prior; their trend target uses the active-donor normalized mean. This
preserves the old reference population rather than mixing heterogeneous residual
DF into its prior-variance estimator. It does not qualify that prior's calibration
for the selected active-donor population.

`supportResolution` records the outcome, original retained/excluded observation
indices, excluded/active donor IDs, rebuilt columns/design/contrast, residual DF,
and active-donor mean. Indices refer to the original contrast observations, not
cells. The original source counts, pseudobulks, global design and size factors
remain in the report. The ordinary feature mean remains descriptive over all
original observations. BH is recomputed across the policy's complete tested
family; original full-support raw inference remains identical, but its adjusted
p-values can change.

## Full real-data results, 2026-09-09

Same archived Kang B-cell and Hagai mouse H5AD matrices, paired designs,
Gamma-trend plans, count filters, offsets and minimum replication as the
[direct Bioconductor comparisons](../../Bioconductor/README.md).

| Dataset | Default tested | Active-policy tested | Insufficient active donors | Remaining support-rank rejection |
| --- | ---: | ---: | ---: | ---: |
| Kang: 2,651 cells, 15,706 genes, 8 donors | 5,400 | 8,829 | 3 | 62 |
| Hagai: 13,863 cells, 22,048 genes, 3 donors | 12,426 | 12,426 | 169 | 78 |

Kang adds 3,429 tested genes with 3–7 active donors. Of its remaining 62 support
rejections, 56 have a zero pair but still lack full positive-support rank after
removal; six have no zero pair. Hagai's zero-pair candidates cannot retain three
donors. Neither dataset has additional numerical-failure or dispersion-boundary
statuses under this option. Original low-expression filters remain 6,812 Kang
and 9,375 Hagai genes.

Both fresh default reports are byte-identical to the previous default reports.
All four fresh product bundles pass source/implementation/replay verification.
Original counts, global designs and prior diagnostics match exactly. All 5,400
Kang and 12,426 Hagai original tested genes retain exact raw fit diagnostics and
unadjusted inference.

Independent NumPy design reconstruction checks every applicable gene against
pseudobulk counts independently aggregated from the archived H5AD. Statsmodels
0.15.0 checks all 3,429 added Kang fits at the reported dispersion: maximum
absolute natural-log effect error 7.48e-6, SE error 5.45e-7, Wald-p error 4.10e-6;
relative mean/influence and absolute leverage errors are below 7.63e-6. Thirty
SciPy adjusted-profile MAP checks (six per active-donor stratum) have maximum
objective gap 8.77e-12. All checks use the predeclared 2e-5 fit tolerance and 1e-4
profile-objective tolerance. No final coefficient-fit reference warnings occurred.
Independent BH checks cover the complete tested families. Hagai has no added fit;
its evidence here is source/design/default preservation and rejection accounting.

### Distinct R methods

Both existing fixed-native-factor and package-normalization edgeR, limma-voom and
DESeq2 tables are compared descriptively. R retains all donor observations and
uses distinct estimators/tests, so this is not a numerical equivalence test.
Across the 3,429 added Kang genes, fixed-factor effect Spearman correlations are
0.97756 (edgeR), 0.98027 (voom) and 0.98584 (DESeq2); sign agreement ranges from
0.92330 to 0.95013. Disagreement remains part of the evidence.

Kang has 815 native BH<0.05 results across 8,829 tests, including 49 added genes.
The default had 850 across 5,400 tests: expanding the tested family changes BH
for the unchanged raw tests. The R comparison JSON separates original published
family values from BH recomputed on each comparison subset. Subset values are
descriptive only and do not replace the native report's family. Among the 65
Kang eligible genes still unavailable natively, fixed-factor R calls 29 (edgeR),
23 (voom) and 6 (DESeq2) significant. Broader coverage and FDR qualification are
not claimed complete.

## Reproduction and evidence

```sh
python run.py --binary /path/to/numivivo \
  --source-bundle /path/to/previous-native-bundle --out /new/path/products
python check.py --baseline /new/path/products/default/report.json \
  --active /new/path/products/active/report.json \
  --counts /path/to/independently-aggregated/counts.tsv --out /new/path/check.json
python compare_r.py --report /new/path/products/active/report.json \
  --reference /path/to/pinned-R-reference --out /new/path/comparison.json
```

The runner explicitly loads native HDF5 from the Apple Homebrew location; adapt
that environment path on other hosts. Python is reference tooling, not part of
the Swift inference runtime. The checker rejects fractional counts; an exact
integer serialized as `1.0` is accepted without rounding.

[evidence/2026-09-09](evidence/2026-09-09) preserves compressed native reports,
plans/receipts, run manifests, per-gene independent checks, R comparisons, source
hashes and logs. Full independent counts and pinned R tables are already retained
in the [Bioconductor evidence](../../Bioconductor/evidence/2026-09-09).
The release binary SHA256 is
`cfceb150c729995af90858c8f2ae382d13bd8b5a0081eb958b7e5c2478396a94`.

All 44 single-cell tests in 12 suites pass on the Mac mini; full release and
scoped H5AD builds and transferred-binary signature/hash checks pass. Analytic fixtures check baseline
rebasing, likelihood/information limits, minimum replication, positive-support
rank and batch confounding; they are numerical regression evidence only.
Preserved unsuccessful checks include the first Swift Testing macro compilation,
an analytic effect assertion tighter than the fitter's 1e-7 scaled-score stopping
criterion, and the reference TSV reader initially requiring integer storage dtype
rather than exact integer values. Corrections and passing reruns are retained;
none changed the solver or rounded source counts to obtain a passing result.
