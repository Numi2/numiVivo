# Native negative-binomial numerical owner

`VivoOmicsNegativeBinomial` implements a Swift NB2 likelihood with variance
`mu + alpha * mu^2`, a log-link GLM with explicit offsets and arbitrary bounded
full-rank covariates, weighted QR Fisher scoring with likelihood step halving,
and Cox–Reid adjusted dispersion profiling. An optional explicit Gaussian prior
on log dispersion supports MAP optimization. No dense cells-by-genes matrix is
used: one gene's donor-level response is fitted at a time.

The cohort CLI still uses the existing log-linear baseline. This numerical
owner is a step toward native NB DE, not a completed replacement. Cohort trend
estimation, empirical dispersion shrinkage/outlier handling, integration with
donor/batch metadata, multiple-testing inference, effect shrinkage and robust
scientific qualification remain.

## Diagnostics and inference boundary

A fit records score convergence, iterations, fitted means, Pearson residuals,
leverage and Cook's distances. Nonconvergence is explicit. Rank-deficient
positive-count support is detected independently of score convergence; effects,
standard errors and adjusted profile likelihood are omitted for these fits.
Dispersion profiling rejects them. This is a conservative support-rank gate,
not a claim to solve every separation or identification case.

Candidate coefficient values and fitted means remain available for numerical
diagnostics; they must not become inferential results when the rank gate fails
or convergence is false. Cook's distances are diagnostic values, without an
outlier threshold, automatic replacement or robustness qualification. The owner
returns no p-values. All coefficient/effect units are **natural logs**.

Dispersion is bounded to [1e-8, 100]. A 25-point log grid locates a bracket,
followed by golden-section refinement. Lower/upper boundary estimates are
explicit. This search is not a proof of global optimality for every arbitrary
input. Inference on raw counts outside the exact FP64 integer range is rejected;
raw exchange retains its separate UInt64 authority.

The adjusted profile method follows the likelihood/information construction
in the [DESeq2 methods documentation](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).
It does not reproduce DESeq2's full dispersion-estimation pipeline.

## Evidence on real counts

[Recorded evidence](evidence/2026-09-09/reference.json) uses the unchanged Kang
pseudobulks from the experimental benchmark: 8,894 eligible genes, 16 donor ×
condition observations and nine design columns. Dispersion **0.15 is supplied
for numerical validation**, not estimated or biologically qualified.

Independent statsmodels GLM fits agree on all 5,400 genes with full-rank
positive-count support. Maximum absolute natural-log effect error was 1.79e-6,
standard-error error 9.97e-8 and log-likelihood error 5.00e-11. The other 3,494
genes are explicitly flagged; they are not counted as successful inferential
fits. NumPy's independent rank check matches every flag. All 293 reference
warnings are retained in the report.

For the first 128 source-order eligible genes, 55 dispersion profiles are
rejected by the support-rank gate, 35 reach the lower bound and 38 interior
profiles agree with an independent SciPy optimizer/statsmodels refit to maximum
objective error 9.80e-8. Near-Poisson boundary optima are **not** claimed to match
that reference: its direct log-gamma subtraction loses precision there.

Fifteen Swift tests passed, including six new tests: analytic intercept mean
and information, offset reparameterization, invalid/nonconverged inputs,
zero-only covariate support, explicit-prior behavior, and eight log masses
calculated independently with mpmath 1.3.0 at 80 digits. Failed development
attempts and their corrections remain in `attempts.json`.

## Reproduction

First reproduce the [Kang benchmark](../Benchmarks/README.md), then use its
output directory as `--kang-result`. The scoped executable calls the actual
Swift numerical owner; it is not a replacement implementation or product CLI.

```sh
python -m pip install -r Tools/Omics/NegativeBinomial/requirements.txt
bash Tools/Omics/NegativeBinomial/build.sh /tmp/numivivo-nb-build
python Tools/Omics/NegativeBinomial/check_reference.py --binary /tmp/numivivo-nb-build/nb-check --kang-result /tmp/kang-result --out /tmp/nb-reference
```

Run Python without `-O`; output directories must be new. The full native input
and output, reference coefficients and warnings are retained locally. The
committed report hashes the executable and exact numerical input. The recorded
native timing is a scoped optimized numerical run, not end-to-end DE throughput
or a speed comparison with Scanpy/PyDESeq2.
