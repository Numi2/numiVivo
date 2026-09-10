# Empirical effect prior for native NB shrinkage

`negativeBinomialOptions.effectPriorEstimation: "weightedUpperQuantile"` learns
a normal-prior width for the requested contrast before the existing count-model
MAP fits. It is mutually exclusive with `effectPriorStandardDeviationLog2`.
Omitting both retains the original unshrunk path; existing fixed-prior plans
keep their meaning. For example:

```json
"negativeBinomialOptions": {
  "trend": "gammaParametric",
  "effectPriorEstimation": "weightedUpperQuantile"
}
```

The estimator matches the weighted 95th percentile of absolute unshrunk log2
effects to the 97.5th percentile of a standard normal. Weights are
`1 / (1 / meanNormalizedCount + trendDispersion)`. Weights are normalized to the
reference-gene count before frequency-quantile interpolation. This follows the
weighted-quantile strategy described in
[DESeq2's prior-estimation implementation](https://github.com/thelovelab/DESeq2/blob/devel/R/core.R).
The independent reference calls the installed DESeq2 1.52.0 routine on the exact
native effects, means and trend dispersions. It does not compare two packages'
different dispersion or coefficient estimates and call them the same prior.

References must be available native Wald fits under the full cohort design,
with absolute log2 effect below 10. At least 20 reference genes are required;
the prior SD has an explicit floor of 0.01 log2 units. Genes with larger effects
are excluded only from prior estimation, not from otherwise eligible shrinkage.
Active-donor fits borrow the full-design prior rather than determining its
width from a mixture of gene-specific donor designs. No p-value threshold is
used to choose reference genes. Existing count, replication, support, dispersion
and optional influence gates still determine availability.

Reports add `negativeBinomial.effectPriorEstimate`, including eligible,
reference and high-effect-excluded feature indices, weight sum, effective
reference count, observed quantile, raw/applied SD and the floor flag.
`effectShrinkage` records the applied width and fit availability. Failure to
estimate a prior produces `effectPriorEstimationError` and no shrinkage summary
or per-gene MAP estimates. It does not silently choose a fixed prior or remove
an original Wald/BH result. The fixed-prior option and the new empirical option
both use the existing joint nuisance-coefficient refit.

The uncertainty remains a conditional Laplace approximation with the estimated
prior treated as fixed. This is a single normal prior, not a heavy-tailed or
mixture prior. Prior uncertainty, effect-estimation risk, posterior coverage,
power and statistical calibration require separate evidence. Learning a width
does not repair the retained Hagai null-calibration concern or the Crowell
reference significance disagreements.

## Reproduction

`prepare.py` freezes ten analyses from the existing fixed-prior Kang, Hagai and
Crowell archives: Kang default and active-donor policies, Hagai default, and all
seven replication-eligible Crowell populations. It retains the full existing
count scopes and records external baseline/source hashes. The original CPE
replication failure remains unavailable. Preparation also writes exact
`prior-input.tsv` tables for the independent R calculation.

```sh
python Tools/Omics/NegativeBinomial/EmpiricalPrior/prepare.py --root RUN --fixed_prior_root FIXED_PRIOR_RUN --crowell_root CROWELL_RUN
Rscript Tools/Omics/NegativeBinomial/EmpiricalPrior/check_prior.R RUN
python Tools/Omics/NegativeBinomial/EmpiricalPrior/run_native.py --root RUN --binary NATIVE_BINARY --hdf5 HDF5_DYLIB
python Tools/Omics/NegativeBinomial/EmpiricalPrior/check.py --root RUN
python Tools/Omics/NegativeBinomial/EmpiricalPrior/archive.py --root RUN --out NEW_EVIDENCE
python Tools/Omics/NegativeBinomial/EmpiricalPrior/archive.py --out NEW_EVIDENCE --verify
```

Use the baseline environment versions in the archived receipts and set
`R_LIBS_USER` for the installed DESeq2 library. Python checks use NumPy/SciPy;
Python and R are reference/preparation tools, not native runtime dependencies.
`check.py --case CASE` checks one completed case; the complete command checks
all ten. No synthetic count data substitute for these experimental runs.

The checker requires exact baseline counts, QC, designs, dispersion diagnostics
and every original feature-inference field. It compares the learned prior with
R, then independently solves all conditional MAP score equations and recomputes
means, likelihoods, penalized objectives and observed-information posterior SDs.
It retains every failed tolerance and solver flag. These are coefficient and
information checks conditional on final dispersion and learned prior, not
evidence of biological truth or calibrated posterior probabilities.

Large complete native reports and original count sources remain externally
retained and hash-bound. The evidence archive contains exact plans, receipts,
reference inputs, per-gene checks and logs, with stored/decompressed checksums.
Use the existing `EffectShrinkage/restore_bundle.py` with the full inventoried
report archive and exact source H5AD for native artifact replay.

## Physical Mac mini qualification, 2026-09-10

The M4 Pro CPU FP64 release build passes all **18 focused Swift tests** across
three suites. The binary SHA256 is
`b0fd73896d17dbfdaadc3173799a56940fe32cd55aa6a2dba852d1b0ad8da59f`.
Ten measured publications and ten successful replays cover **99,394 MAP fits**.
Every conditional numerical check passes; no SciPy solver failure flags occur.
All original counts, QC, designs, dispersion diagnostics and unshrunk feature
inference remain exactly equal to the archived baselines.

| Analysis | MAP fits | Prior reference genes | Learned SD, log2 |
| --- | ---: | ---: | ---: |
| Kang default | 5,400 | 5,400 | 1.126058 |
| Kang active donors | 8,829 | 5,400 | 1.126058 |
| Hagai default | 12,426 | 12,412 | 1.450802 |
| Crowell astrocytes | 10,650 | 10,650 | 0.683190 |
| Crowell endothelial | 10,140 | 10,140 | 1.085796 |
| Crowell microglia | 9,435 | 9,435 | 0.917023 |
| Crowell oligodendrocytes | 10,387 | 10,387 | 0.537446 |
| Crowell OPC | 10,074 | 10,074 | 0.562134 |
| Crowell excitatory neurons | 11,070 | 11,070 | 0.286484 |
| Crowell inhibitory neurons | 10,983 | 10,983 | 0.264780 |

The largest native/R prior-SD difference is **4.22e-15 log2 units**. The largest
native/SciPy natural-log MAP effect difference is **7.54e-8**, posterior-SD
difference **7.84e-9**, relative mean difference **2.21e-7**, likelihood difference
**2.56e-7** and objective difference **6.05e-9**. Every recomputed scaled score
is below the unchanged 1e-7 tolerance. Fourteen large-effect Hagai genes are
excluded from prior estimation while retaining eligible MAP fits. All 3,429
additional active-donor Kang fits use the full-design prior, with their own
retained observations and nuisance-coefficient refits independently checked.

One excitatory-neuron replay attempt failed with `omics snapshot write` under
low disk space. The failed log and command remain retained; after removing
verified regenerable compiler cache files, the same published bundle replayed
successfully with the unchanged binary. The driver resumes from checked
successful artifacts and retains all attempts rather than restarting the suite.
This is 21 native command attempts for the ten empirical-prior cases, including
that one failure. Source data, binaries and research evidence were preserved.

A separate compatibility check republishes the unchanged fixed-prior Kang
plan and successfully replays it with the new binary. Its entire plan and
report are **byte-identical** to the previous fixed-prior outputs. Direct replay
of the old receipt with the new implementation correctly fails its fingerprint
guard; that rejection is retained and the receipt is not rewritten. Numerical
compatibility and immutable implementation-bound receipt validity are separate.

[Evidence](evidence/2026-09-10/manifest.json) retains complete per-gene checks,
all native attempt logs, source/protocol hashes and R versions. These are ten
analyses of three studies with overlapping genes and observations, not ten
independent studies. Native build logs retain pre-existing warnings from
unrelated Metal/workflow code. Neither this prior estimator nor its conditional
fit agreement qualifies FDR, posterior coverage or a general accuracy gain.
