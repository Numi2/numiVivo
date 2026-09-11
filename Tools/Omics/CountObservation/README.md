# Conditional uncertainty from sampled RNA counts

`VivoCountObservation` adds the count observation layer needed to address the
[GSE181897 sampling diagnosis](../PerturbationPrediction/GSE181897/Uncertainty/README.md).
It evaluates a single gene's latent population rate from original cell counts
and their complete RNA library sizes. Zero observations retain posterior
uncertainty under an explicit proper prior. It also computes moments for a
planned new sample of the same condition, separating measurement variation from
uncertainty in that latent rate.

This is a native statistical component. It does **not** repair or replace the
published perturbation intervals, infer dispersion, learn a biological prior,
or establish calibrated prediction of a treated donor. Those connections and
independent biological validation remain open.

## Count model and units

For cell `i`, full-RNA library size `L_i`, gene count `Y_i`, and latent rate
`r` in CPM:

- `E[Y_i | r] = L_i * r / 1e6`.
- `Var[Y_i | r] = E[Y_i | r] + phi * E[Y_i | r]^2`.
- `r ~ Gamma(shape, ratePerCPM)`, with both prior parameters supplied explicitly.
- Conditional cells are independent. `phi=0` is the exact Poisson limit.

NumiVivo's `phi` is dispersion; it is the reciprocal of the precision parameter
in the [Stan NB2 definition](https://mc-stan.org/docs/functions-reference/unbounded_discrete_distributions.html#negative-binomial-distribution-alternative-parameterization).
The [Gamma distribution](https://mc-stan.org/docs/functions-reference/positive_continuous_distributions.html#gamma-distribution)
is parameterized by shape and rate, not scale. The supplied prior is part of
the model, not an epsilon added to zero counts. There is no default prior or
dispersion. A prior dominated result is not data-driven evidence of expression.

The posterior includes means and variances in CPM and natural-log(1+CPM), plus
equal-tail latent-rate credible bounds in natural-log(1+CPM). These bounds
are **not** future observed-count intervals. Rates are marginal NB regression
rates; their support is unbounded and genes are not jointly constrained to a
closed RNA composition. Library sizes are conditioned-on offsets, not modeled
random denominators. Correlated cells, cell-state mixtures, donor variation,
assay differences and uncertainty in fitted dispersion are not resolved merely
by supplying a positive dispersion.

Each count/depth pair must represent an actual independent cell at this API
boundary. A pre-summed donor row is not one cell, and a source-cell count alone
cannot reconstruct the likelihood at unequal cell depths. Sparse callers can
extract one gene at a time while keeping the full-RNA denominator. The statistical
function does not inspect source identities: the caller must bind cell membership,
feature/condition identity, source counts, prior and dispersion provenance. The
benchmark below records those bindings explicitly.

## Planned sampling

`predictiveMoments` requires each planned cell's complete RNA depth. With
`e_i = plannedL_i / 1e6`, `E = sum(e_i)`, posterior rate mean `m` and variance `v`:

- Mean gene counts: `E*m`.
- Conditional Poisson variance: `E*m`.
- Conditional cell overdispersion: `phi * sum(e_i^2) * (v + m^2)`.
- Latent-rate contribution: `E^2*v`.

Their sum is the future count variance by total variance. At equal total depth,
more independent cells reduce the cell-overdispersion term. The latent-rate term
is shared across future cells and must not be divided by the planned cell count.
This endpoint is a new independent sample of the **same** latent condition;
no perturbation response is applied. These moments do not supply a count
predictive interval or moments of log-normalized future counts.

## Numerical method

A proper Gamma prior makes the posterior in log-rate strictly log-concave.
The implementation brackets its unique mode, uses concavity-based tangent bounds
to bound omitted tails, and integrates with composite 16-point Gauss-Legendre
quadrature. It doubles panels until posterior moments and quantiles stabilize.
It integrates centered moments to retain small variances at high counts and
uses Taylor remainders to avoid cancellation around the mode. Invalid input,
failed tail bounds and nonconvergence produce errors, never successful intervals.

The posterior reports the refinement difference, number of panels and an
analytic tail bound for probability and scaled first/second rate moments.
The refinement check is numerical evidence, not a rigorous discretization error
bound or biological coverage guarantee. The output does not report the marginal
likelihood or posterior density normalizer.

## Qualification and reproduction

The numerical protocol uses all 62 admitted GSE181897 **control** donors,
2,309 original cells and 16 panel genes selected by SHA256 order of their IDs.
Dispersion values 0, 1 and 10 are fixed sensitivity cases, with an explicitly
unqualified Gamma(0.5, 0.001 per CPM) prior. The resulting 2,976 cases include
1,503 all-zero gene/donor cases. Planned sampling is 20 cells at 1,000 RNA UMIs
each, fixed without reading treated count values or treated depths. All 62
control aggregates are compared with the previously qualified full-RNA counts.

These are numerical checks using real counts, **not** 2,976 independent biological
replicates or a fitted sampling model. GSE181897 is already known development
data. Seventy-five additional numerical boundary cases cover zeros, depth
imbalance, high counts, near-Poisson dispersion and broad proper priors.

The reference uses independent adaptive QUADPACK integration, with 60-digit
likelihood evaluation for boundary cases. Poisson cases also have an exact
Gamma-conjugate posterior check. Native output is evaluated for posterior
moments, credible bounds and planned-sample variance components.

```sh
bash Tools/Omics/H5AD/build.sh /absolute/new/build --with-cli
bash Tools/Omics/CountObservation/test.sh /absolute/new/build
python Tools/Omics/CountObservation/prepare.py /path/to/original/GSE181897 /new/inputs
/absolute/new/build/count-observation-check /new/inputs/cases.json /new/controls.jsonl
python Tools/Omics/CountObservation/verify.py /new/inputs/cases.json /new/controls.jsonl /new/reference.json
```

The preparation requires AnnData, NumPy and SciPy and the original parent study.
Verification needs NumPy, SciPy and mpmath but can run from the retained exact
cell vectors without the original H5AD. `Check.swift` is a qualification runner
calling the product API, not an added production CLI command.

## Completed results

On an Apple M4 CPU with Swift 6.3, all **25 tests in three suites** pass,
including the existing NB2 and perturbation regression tests. All **2,976
experimental-count cases and 75 numerical boundary cases** pass the independent
reference checks at the predeclared 2e-7 tolerance. The largest relative moment
error is 7.10e-8. All 1,503 zero-count cases have positive posterior variance.
The full control run repeats byte-for-byte.

Native control execution takes 4.37 seconds with 31,113,216 bytes (29.67 MiB)
peak RSS. This is a measured qualification run, not a controlled speed comparison
against SciPy. Native panels are at most 64 for these control cases, and the
largest reported scaled-moment tail bound is 5.02e-16. Independent reference
checks take 48.44 seconds for controls and 12.38 seconds for boundary cases.

The [results](evidence/2026-09-11/results.json) and
[complete archive manifest](evidence/2026-09-11/manifest.json) retain every
original cell vector, both native control outputs, boundary output, independent
reference result, protocol, source hashes, recipes, environment and execution
logs. The archive can be extracted with `tar -xzf results.tar.gz`; no original
H5AD is needed to rerun the retained numerical cases. The original H5AD remains
required to regenerate and verify those cell vectors from the source.

An initial test-harness type error and two false convergence failures are
retained. The latter exposed a comparison of an analytically zero centered
moment correction; the final method compares the actual posterior mean.
The rebuilt implementation passes all checks. No full package test, GPU
acceleration, learned-prior assessment or biological coverage claim is implied.

## Next integration work

Fit or specify source-bound cell dispersion and rate priors using admitted
training/control data, preserving uncertainty in those estimates. Then connect
observation likelihoods to a latent donor response model, separating training
measurement noise from donor variation rather than adding it twice. Query
prediction must require the intended latent or future sampled endpoint and,
for the latter, planned sampling information. Preserve the original failed
interval scores and validate the combined model on new biological observations.
