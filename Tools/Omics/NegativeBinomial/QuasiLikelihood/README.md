# Quasi-likelihood stages and native deviance

This study resolves the inputs and arithmetic needed for native QL development.
It does not introduce a native QL cohort test. Current Wald and LRT reports and
options remain unchanged. The new production primitive is
`VivoOmicsNegativeBinomial.unitDeviance(count:mean:dispersion:)`, the residual
likelihood calculation needed before QL scale estimation and moderation.

The frozen [protocol](PROTOCOL.md) uses all twenty original default-support
Kang/Hagai sham analyses and nine full-support treatment analyses, including
seven Crowell populations. All source counts, feature identities, offsets,
designs and stored native trends pass an independent SciPy CSR re-read of the
hashed source reports. Active-donor profiles are explicitly outside this study:
borrowing a QL prior across their differing designs needs its own contract.

Each case runs four edgeR 4.10.5 fits with robust QL prior estimation: native
stored trend or edgeR estimated trend, each with modern adjusted or legacy
QL. None uses native gene-wise MAP dispersions as a QL trend. This follows the
separation of NB trend and residual quasi-dispersion described in the
[edgeR reference manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf).
Offsets differ from native relative factors by one common effective-library
constant, with the same offsets in all arms. No TMM or new filtering is applied.
The inspected package functions, arguments and recorded solver diagnostics
identify the actual reference path; no native implementation is inferred from R.

## Results on 2026-09-10

All 116 reference fits and exact-dispersion method constraints pass. Across
967,152 gene/analysis/arm observations, independent F, denominator DF, posterior
scale, F-tail/Poisson-bound and BH checks pass. Adjusted unit deviances and DF
sum to their published aggregates. Maximum p-value error is 1.28e-15 and BH
error 3.56e-15. No package failure flags occur, but maximum scaled score is
5.61e-4; a false failure flag is not a strict convergence certificate.

| Study | Native trend adjusted QL | Native trend legacy QL | edgeR trend adjusted QL | edgeR trend legacy QL |
| --- | ---: | ---: | ---: | ---: |
| Kang sham BH<0.05 calls | 2 (1/10 splits) | 0 (0/10) | 3 (2/10) | 0 (0/10) |
| Hagai sham BH<0.05 calls | 39 (8/10 splits) | 38 (8/10) | 79 (10/10) | 83 (9/10) |

These are gene/split call totals on overlapping, previously inspected sham
splits, not independent experiments. The same native default families previously
had Wald/LRT totals of 0/0 for Kang and 20/22 for Hagai. QL does not establish a
calibration repair here. Neither fewer calls nor package agreement selects a
production default. Earlier package benchmarks admitted additional genes; their
call totals are not interchangeable with this common-family experiment.

| Treatment | Genes | Native trend adjusted / legacy | edgeR trend adjusted / legacy |
| --- | ---: | ---: | ---: |
| Kang | 5,400 | 906 / 931 | 1,072 / 1,046 |
| Hagai | 12,426 | 4,613 / 4,562 | 4,904 / 4,812 |
| Crowell Astro | 10,650 | 466 / 477 | 497 / 509 |
| Crowell Endothelial | 10,140 | 569 / 574 | 625 / 620 |
| Crowell Microglia | 9,435 | 219 / 215 | 281 / 274 |
| Crowell Oligodendrocytes | 10,387 | 56 / 64 | 65 / 79 |
| Crowell OPC | 10,074 | 124 / 134 | 153 / 153 |
| Crowell Excitatory neurons | 11,070 | 32 / 13 | 33 / 13 |
| Crowell Inhibitory neurons | 10,983 | 6 / 5 | 5 / 5 |

Treatment counts are BH<0.05 discoveries, not known power or effect truth.
All per-gene results and diagnostic distributions are retained.

## Native arithmetic

The native primitive evaluates twice the saturated-minus-model log likelihood.
It supports exact integer counts through 2^53, nonnegative finite means, zero
count/zero mean, the exact Poisson limit, and NB dispersions in [1e-8,100].
Positive counts at zero mean, nonfinite inputs, unsupported dispersion and
unrepresentable outputs fail explicitly. Near saturation it evaluates the
entropy difference as one convergent series, keeping the small dispersion
factor inside each term; it does not subtract two nearly equal log masses.
Other branches use scaled log ratios and explicit zero-count formulas.

On the physical M4 Pro CPU, all 22 focused Swift tests pass. All 7,881,944 actual
reference fitted-count pairs agree with independent direct likelihood arithmetic
within 6.36e-10 absolute error, without high-precision refinement. The fixed
702-case boundary grid agrees with 100-digit calculations within 7.48e-14
relative error, including tiny representable deviances at counts near 2^53.
Summed native deviances differ from edgeR's stored raw deviances by up to
2.64e-6; this additional package-arithmetic observation is retained separately
from the independent calculation at identical means.

The qualified native executable SHA256 is
`db0744334d031d76bff070d1b1e0e9bd2dedd9cf43b0c35dda58cb2a257ee72f`.
This is CPU FP64 primitive validation, with no Metal or performance claim.

## Next native ownership

[Native direct moments and adjusted residuals](Moments/README.md) now implement
conditional unit-deviance moments and leverage-dependent residual DF at supplied
means, trend dispersion and average QL scale. All 3,940,848 available default-budget
moments pass independent checks. Twenty gene/arm fits exhaust the default work
limit; a separate higher-work sensitivity resolves nineteen, leaving one Hagai
fit unavailable. The subsequent [explicit adaptive evaluator](Moments/Adaptive/README.md)
closes all twenty measured work-limit failures: all 58 arms and 3,940,972 moments
pass, including independent checks of the 124 recovered moments. Direct remains
the default, and the original failures remain retained.

A native QL cohort method still needs the global QL scale/refit, robust abundance-
dependent prior estimation for unequal residual DF, and constrained-test
integration with explicit failure behavior. The existing equal-DF, untrended,
non-robust linear variance prior is insufficient. The conditional residual stage
must not be exposed as if it provided those missing stages. The reference
matrices here supply stage-level comparison targets, including unit adjustments,
prior/posterior scales and optimizer diagnostics. Scientific calibration,
independent power, uncertain-dispersion coverage and varying-support designs
remain separate requirements after implementation.

## Reproduction and evidence

```sh
python prepare.py --lrt-root PRIOR_LRT_RUN --out RUN
python run.py --root RUN --r-library PINNED_R_LIBRARY --jobs 2
python check_inputs.py --root RUN
python check.py --root RUN
bash build.sh NATIVE_RUN/build
python check_native.py --root RUN --binary NATIVE_RUN/build/nb-deviance
python archive.py --root RUN --out ARCHIVE
python archive.py --out ARCHIVE --verify
```

Reference orchestration uses pinned edgeR 4.10.5, limma 3.68.5, statmod 1.5.2
and jsonlite 2.0.0. Independent checks use NumPy/SciPy and mpmath; native builds
use Swift 6.3.3. `check_native.py` sends donor-sized JSON batches to the actual
production-owner executable on the source manifest's host. No dense cell-by-
gene matrix is constructed.

The [archive manifest](evidence/2026-09-10/manifest.json) binds all pseudobulk
inputs, all four complete QL result tables per case, source and arithmetic
checks, native check summaries, boundary outputs and test/build logs. Complete
reference fit matrices and native per-observation outputs remain externally
retained at the listed paths and hashes. They are not all stored in Git. The
[reference summary](evidence/2026-09-10/summary.json) and
[native summary](evidence/2026-09-10/native-summary.json) distinguish the two
owners and their evidence scopes.
