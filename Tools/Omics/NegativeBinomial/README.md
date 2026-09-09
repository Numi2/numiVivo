# Native negative-binomial cohort analysis

The actual `singlecell-analyze` route now supports an explicit
`"model": "negativeBinomial"` contrast. It uses the existing sample selection,
paired-donor/independent-replicate checks, batch design, sparse pseudobulk and
size-factor authority. Old plans without a model retain the log-linear baseline.
The NB method remains experimental pending multi-study calibration and broader
reference coverage; enabling it does not claim production qualification.

## Model and diagnostics

`VivoOmicsNegativeBinomial` implements NB2 with variance
`mu + alpha * mu^2`, explicit log size-factor offsets and a log-link GLM.
Observed-information Newton steps use weighted QR and likelihood step halving;
final covariance uses expected Fisher information. Coefficient convergence,
fitted means, Pearson residuals, leverage and available Cook's distances are
retained. Core coefficients use natural logs; cohort effects, standard errors
and intervals use log2 units.

Gene-wise dispersion maximizes the Cox–Reid adjusted profile likelihood.
A robust Huber fit of log-dispersion residuals estimates the parametric trend
`alpha = a0 + a1 / mean`, with nonnegative coefficients. The explicit `mean`
alternative uses the median log dispersion and never activates as a silent
fallback. Only interior gene-wise estimates enter trend fitting.

The explicit `gammaParametric` alternative fits the same mean curve using
the Gamma identity-link objective `sum(alphaGene / fitted + log(fitted))`.
Positive coefficients are fitted by Fisher scoring with a checked line search.
After each fit, genes with observed/fitted dispersion below 1e-4 or at least
15 are excluded from the next trend fit. Identifiability, positive coefficients,
minimum retained genes, inner stationarity and outer convergence are required;
failure does not switch methods. `trendFitFeatureIndices` records the retained
subset, while `referenceFeatureIndices` still records the full interior cohort
used for the prior variance calculation. Existing `parametric` and `mean`
plans retain their previous semantics and omit the new optional field.

This separates a Gamma mean trend from the robust trend of log dispersions.
It does not reproduce the entire DESeq2 procedure: native gene-wise dispersion
refits coefficients along the profile, whereas the tested PyDESeq2 version
uses its preliminary fitted means during dispersion estimation. Prior estimation
and feature admission also remain explicit native choices.

The log-normal prior variance is the squared scaled MAD of trend residuals
minus `trigamma(residualDF / 2)`, floored by the declared minimum. MAP estimation
then refits each gene's adjusted profile with this prior. High positive
log-dispersion residuals beyond the declared outlier threshold retain their
gene-wise estimate. The trend, reference genes, prior, outliers, original/final
dispersions and per-gene failures are all present in the report.

Normal Wald probabilities and intervals condition on the final dispersion;
they do not integrate all uncertainty in trend/prior estimation. BH adjustment
uses only the available tests within a contrast. Repeated-sampling calibration,
selection-adjusted FDR and cross-contrast guarantees are not established.
Very small probabilities can underflow to zero in FP64; zero is not a claim
of mathematical impossibility. The recorded native and PyDESeq2 p-values
sometimes differ by many orders of magnitude.

Positive-count support rank deficiency withholds inference, even if candidate
coefficient optimization reports score convergence. This conservative gate
currently excludes 3,494 Kang genes and is not a complete treatment of every
gene-specific nuisance-parameter or separation case. Failed fits, final boundary
estimates and genes excluded by an explicitly requested influence threshold
receive no p-value. Cook's distances are unavailable for deficient support or
unit leverage; they are never replaced with fabricated zeros. No count
replacement or automatic donor removal occurs.

Dispersion search is bounded to [1e-8, 100], using a 25-point log grid and
refinement. Estimates within 1e-4 log units of an endpoint are conservatively
flagged as boundary estimates, without replacing their recorded values. This
accounts for near-Poisson likelihood rounding beyond the scalar search width.
The search is not a global-optimality proof for arbitrary inputs.

The adjusted-profile construction is described in the
[DESeq2 methods documentation](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).
Our robust trend, optimization, rank gate and diagnostics are explicitly
reported choices; this is not a reproduction of the entire DESeq2 pipeline.

## Native plan

Add these fields to an otherwise complete existing contrast:

```json
{
  "model": "negativeBinomial",
  "negativeBinomialOptions": {
    "trend": "parametric",
    "minimumTrendGenes": 20,
    "minimumPriorVariance": 0.25,
    "outlierStandardDeviations": 2
  }
}
```

An optional positive `maximumCooksDistance` excludes genes above that declared
threshold. Omission reports influence without exclusion. No universal threshold
or robustness guarantee is asserted. Overrides of the log-linear `variance`
and `priorCount` options are rejected for NB rather than silently applied.

`singlecell-analysis-export` reconstructs the complete NB result from its
count/plan receipts. `singlecell-analysis-tables` exports an NB diagnostics JSON
and a `z` column, leaving `t` and `df` empty. Its index identifies the actual
method. Filtered/unavailable probabilities remain empty, not zero.

## Real-data evidence

The [cohort evidence](evidence/2026-09-09-cohort) records a full optimized CLI
run on all 2,651 Kang B cells, all 15,706 source genes and eight paired donors.
Counts, pseudobulks and Scanpy QC match exactly; normalization error is below
1e-10. Of 8,894 expression-eligible genes, 5,400 are tested and 3,494 fail the
positive-support rank gate. There are no numerical failures. The parametric
trend uses 2,823 interior genes; 37 dispersion outliers retain gene-wise fits.

On the 5,400 common tested genes, native versus PyDESeq2 effect Spearman
correlation is 0.99907, sign agreement is 99%, and 43 of the top 50 BH-ranked
genes overlap. All five declared interferon genes pass the direction check.
These are descriptive results on one selected cell type in one study, not
competitiveness or calibration gates. PyDESeq2's mean-trend fallback warning
and disabled Cook's/independent-filtering policy remain explicit.

Independent SciPy/statsmodels checks cover the actual fitted trend, prior,
all 5,400 available Wald tests, 32 MAP profiles and BH calculation. Maximum
log2 effect difference is 1.12e-5 and MAP objective difference is 3.15e-12.
The complete CLI also passed NB table/diagnostic checks, 18 Swift tests and
13 baseline CLI assertions across 20 commands, including workflow replay.
Failed development attempts remain in `attempts.json`.

The [earlier conditional-owner evidence](evidence/2026-09-09) and fresh
`newton-reference.json` isolate numerical fitting at supplied dispersion 0.15;
that supplied value is not a biological estimate. Near-Poisson reference
profiles remain explicitly unqualified because of reference log-gamma precision.

## Reproduction

Use Python 3.12 and the pinned requirements. Fresh output directories are
required; run Python without `-O` so integrity assertions remain enabled.

```sh
python -m pip install -r Tools/Omics/NegativeBinomial/requirements.txt
python Tools/Omics/Benchmarks/run_kang.py --binary /path/to/numivivo --source kang_2018.h5ad --model negativeBinomial --nb-trend parametric --out /tmp/kang-nb
python Tools/Omics/NegativeBinomial/check_cohort.py --kang-result /tmp/kang-nb --out /tmp/kang-nb-reference
/path/to/numivivo singlecell-analysis-tables /tmp/kang-nb/analysis-receipt.json --store /tmp/kang-nb/store --output /tmp/kang-nb-tables
python Tools/Omics/NegativeBinomial/check_tables.py --report /tmp/kang-nb/native-report.json --tables /tmp/kang-nb-tables --out /tmp/kang-nb-tables-check.json
```

For isolated numerical development, `build.sh` creates `nb-check` from the actual
owners. Its `cohort dataset.json plan.json report.json` mode is a scoped harness,
not the product CLI. No Python implementation supplies native fitting.

## Remaining work

Multi-study calibration, robust reference sensitivity, gene-specific nuisance
handling and effect shrinkage remain. [Direct Bioconductor comparisons](../Bioconductor/README.md)
now cover Kang and Hagai with fixed native and independent package normalization.
Native support-rank coverage and calibrated significance remain open.
Do not promote the default on effect correlation alone. The original Hagai
32.85-million-nonzero scope now passes through the streaming count route and a
[three-pair LPS6 benchmark](../Benchmarks/README.md#hagai-paired-negative-binomial-comparison).
All 12,426 tested native Wald fits, the trend/prior, 32 MAP objectives and BH
arithmetic passed independent checks. Native/reference final dispersions have a
median ratio of 0.37368, despite effect correlation 0.99919. The resulting
significance differences and the two residual degrees of freedom make this
direct evidence for further calibration work, not permission to promote the
experimental model. Source-prefix donor pairing is explicitly documented.

`check_cohort.py --cohort-report ... --counts ... --out ...` accepts a streaming
report with one source sample per pseudobulk and checks sample ordering before
independent fitting; `--kang-result` remains supported for the resident route.
`compare_dispersions.py` records estimator differences and reference convergence
flags without excluding inconvenient genes or retuning either estimator.

## Explicit Gamma trend qualification

Use `negativeBinomialOptions: {"trend": "gammaParametric"}` in a contrast, or
`--nb-trend gammaParametric` in the Kang/Hagai benchmark scripts. The default
remains the previously recorded robust log trend. The
[Gamma evidence](evidence/2026-09-09-gamma/) holds both original and candidate
comparisons rather than replacing the earlier results.

| Complete experimental scope | Robust log trend | Gamma trend |
| --- | --- | --- |
| Kang, 5,400 tested genes | Effect correlation 0.999067; top-50 overlap 43; 964 BH values below 0.05 | Correlation 0.999331; overlap 44; 850 BH values below 0.05 |
| Hagai, 12,426 tested genes | Effect correlation 0.999190; top-50 overlap 39; 6,283 BH values below 0.05 | Correlation 0.999869; overlap 45; 4,827 BH values below 0.05 |

Correlations and top-50 overlaps compare each native choice to the same
PyDESeq2 reference on commonly tested genes. Significance counts are descriptive
within each native hypothesis family, not false-discovery calibration. Counts,
QC, metadata, experimental design and every contrast option except the trend
are identical across native comparisons. Both Gamma runs retain all five
preselected positive response genes and have zero numerical failures.

Hagai's Gamma trend retains 9,312 of 9,377 interior estimates for curve fitting.
Its intercept is 0.0332723 and inverse-mean coefficient 2.90502. The median
native/PyDESeq2 final-dispersion ratio increases from 0.37368 to 0.72605;
remaining disagreement is recorded. Nineteen Swift tests and the existing CLI
suite pass. Independent statsmodels/SciPy checks cover all 17,826 tested Wald
fits, both trend objectives and priors, 64 MAP objectives and both BH results.
Full-product report hashes match those independently checked reports, and both
bundles reconstruct. This adds an explicit estimator choice with experimental
and numerical evidence; calibrated FDR and production promotion remain open.
