# Native contrast likelihood-ratio test

Set `negativeBinomialOptions.testMethod` to `"likelihoodRatio"` in an existing
negative-binomial contrast plan to test the requested contrast against zero.
Omitting the option retains the original Wald inference and report encoding.
The count/analysis/H5AD publication and replay routes share this production
cohort implementation.

The test preserves the full-model dispersion estimate, raw counts, offsets,
donor/batch design and support policy. It refits nuisance coefficients under
the single linear constraint, compares unpenalized NB likelihoods and uses the
asymptotic chi-square tail with one degree of freedom. The null model is built
from the contrast, rather than assuming treatment occupies a particular column.
Count-only likelihood constants are cancelled before computing the ratio.
This follows the fixed-dispersion GLM likelihood-ratio comparison described in
the [edgeR user guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).

Feature `pValue` and `adjustedPValue` now represent LRT and its own available
BH family when selected. `zStatistic`, standard errors and confidence intervals
remain Wald diagnostics; they are not inversions of the likelihood-ratio test.
Reports identify the method explicitly and retain the null fit under
`negativeBinomial.features[].likelihoodRatioFit`: original-coordinate null
coefficients, means, log likelihood, convergence, scaled score, signed raw LR,
one degree of freedom, statistic and probability. Null optimization failures
withhold inference and retain diagnostic errors. A signed LR below -1e-7 is
unavailable; smaller negative roundoff is recorded and clamped to zero.

Existing count, donor replication, dispersion boundary, support and optional
influence gates remain. Active-donor requests use the same gene-specific rows
and columns for both full and constrained fits. An effect prior does not enter
either likelihood. Optional effect shrinkage remains a separate calculation.

This is an experimental alternative, not a repair established by the earlier
null benchmark. Both Wald and this LRT condition on estimated dispersion;
neither accounts for its uncertainty through quasi-likelihood. Previously
inspected treatment data do not provide effect truth or power, and overlapping
sham splits do not provide independent experiments or general FDR calibration.
No default is promoted by a lower sham-call count or numerical agreement.

## Reproduction

The frozen [protocol](PROTOCOL.md) binds 40 original Kang/Hagai null analyses
and all ten prior treatment analyses across Kang, Hagai and Crowell. Python
orchestrates and independently checks native production-owner calculations;
the R adapter uses pinned edgeR 4.10.5 with native dispersions and offsets and
`prior.count=0`, rather than estimating a different dispersion model.

```sh
python prepare.py --out RUN
bash build.sh NATIVE_RUN/build
python run_native.py --root RUN --binary NATIVE_RUN/build/nb-lrt --jobs 4
python check.py --root RUN --refine-reference
python run_reference.py --root RUN --r-library PINNED_R_LIBRARY
python run_reference.py --root RUN --r-library PINNED_R_LIBRARY --tighter-full-fit
python summarize.py --root RUN
python archive.py --root RUN --out ARCHIVE
python archive.py --out ARCHIVE --verify
```

Paths in the frozen source manifest refer to retained native reports on the
owning Mac mini. The runner verifies source and executable hashes and streams
compact output, avoiding copies of full source data or dense cell matrices.
`check.py` verifies compressed and decoded source hashes before using original
pseudobulk counts. The native harness recomputes original requests before the
new method, reporting exact preservation checks and every full-model fit.

## Measured result: 2026-09-10

On the physical M4 Pro Mac mini, all 21 focused Swift tests and all 50 native
analyses passed. Independent checks cover 443,705 tested gene/analysis pairs:
the original requests, full-model fits, dispersions and available families
remain exact, with no unavailable null tests. The maximum LR arithmetic error
is 7.92e-9, p-value error 9.87e-9 and BH error zero. Two real H5AD publications
and two native replays passed. The unchanged Wald report is byte-identical to
its archived baseline; the product LRT features and null fits exactly match
the measurement harness.

| Sham study and policy | Wald BH<0.05 calls | LRT BH<0.05 calls | LRT splits with calls |
| --- | ---: | ---: | ---: |
| Kang default | 0 | 0 | 0/10 |
| Kang active donors | 0 | 1 | 1/10 |
| Hagai default | 20 | 22 | 8/10 |
| Hagai active donors | 20 | 22 | 8/10 |

These are totals of gene/split calls, not distinct genes or independent
experiments. The alternative does not resolve the sham-discovery concern.
The ten treatment comparisons cover 99,394 gene/analysis pairs; changed
discovery counts do not establish increased power or biological truth.

Default edgeR comparisons pass for all seven Crowell populations but retain
35 p-value disagreements across the three Kang/Hagai analyses (maximum
6.67e-5 versus the declared 2e-6 tolerance). The separately declared
[stopping sensitivity](SENSITIVITY.md), applied through the actual
[solver adapter](SENSITIVITY_ADAPTER.md), passes all ten cases: maximum LR
difference 6.02e-6 and p-value difference 3.33e-8. No reference BH<0.05 decision
changes. This supplementary result does not overwrite the original failures.
The first attempted sensitivity passed options that pinned `glmFit.default`
ignores; all ten outputs remain recorded as a failed method constraint.

The initial independent checker also retained two failed sham checks caused
by cancellation in this platform's 64-bit `longdouble`. An 80-digit calculation
at the exact saved means confirmed native accuracy. The final checker refines
nine candidates at higher precision, without changing native arithmetic or
acceptance tolerances. Original checks, build/storage failures, the repaired
null-helper reporting-contrast failure and an initial missing-HDF5-path attempt
are retained alongside successful results.

The [archive manifest](evidence/2026-09-10/manifest.json) binds final per-gene
checks, reference tables, original check summaries, product reports/receipts
and failed attempts. Complete native measurements, initial per-gene tables,
source model archives and qualified executables are retained externally at
the exact paths and hashes listed in the manifest; they are not all embedded
in Git. The [summary](evidence/2026-09-10/summary.json) deliberately reports
`completed-with-retained-failures`. This is CPU FP64 numerical and bounded
product validation, without Metal, general calibration or production claims.
