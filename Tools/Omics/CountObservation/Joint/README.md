# Joint count-response likelihood: development qualification

**Finite-grid likelihood and conditional prediction arithmetic pass. Grid
refinement fails for 17 of 19 available origin/gene models. Biological accuracy
and uncertainty are not qualified by this experiment.**

This implements the next step after the [paired-donor diagnosis](../Paired/README.md):
fit a coherent nonnegative joint distribution of control and treated RNA rates
from paired training donors, then update it using only a new donor's control
counts. The implementation is `VivoJointCountResponse`, using
`VivoCountRateLikelihood` and the existing count-sampling variance decomposition.
It is a research API and scoped runner, not a qualified production CLI.

## Inputs, endpoint and separation

The original Kang and HIRISA sparse H5ADs supply all 122,164 previously admitted
training cells: 2,651 cells from eight Kang donors and 119,513 from five HIRISA
donors. Every full RNA count contributes to each cell's denominator. All 26
complete group count vectors match the published [calibration](../Calibration/README.md).
The model is executed only for the same 16 identifier-hash-selected development
genes as the previous query experiment. This is not full-transcriptome model
execution or a newly selected successful subset. Unavailable genes remain in
all inputs, outputs and counts.

For each condition, the endpoint is a latent rate `r` in CPM with
`E[Y_i | r] = L_i * r / 1e6`, where `L_i` is the cell's full RNA count. Conditional
counts follow NB2 with the previously estimated fixed cell dispersion `phi`.
A constant rate within a donor/condition is a modeling assumption; when it is
violated, this endpoint need not equal either the arithmetic mean of per-cell
CPM or RNA-weighted pseudobulk CPM.

The query package contains original control cells from 62 GSE181897 donors and
a declared future sample of 20 treated cells with 1,000 RNA counts each. No query
treated outcome is used or scored. The native fit function receives only paired
training strata; query values do not choose support or weights. The runner's
source fingerprint binds the complete frozen origin package, including its
query section, and should not be misread as separately hashed training and query
packages. This is development on known studies, not new biological validation.

## Exact grouping and finite-support estimation

Only cells with exactly equal full RNA depth are grouped. At fixed dispersion,
the rate-dependent log likelihood is

```
Y_total * log(r)
  - sum_depth (Y_depth + N_depth / phi) * log(1 + phi * depth * r / 1e6).
```

The omitted terms do not depend on rate; they can depend on dispersion. This
reduction is sufficient for rate fitting at fixed phi, not for re-estimating
phi. The Poisson limit is explicit. Zero-rate support is allowed for all-zero
counts, and positive counts at zero rate have zero likelihood. The implementation
uses centered stable relative likelihoods and explicit numeric UInt64 conversion.

The two rate axes span the training donor maximum-likelihood rates, with nested
log1p-spaced grids of 9, 17, 33 and 65 points plus every donor MLE. Their Cartesian
product forms a declared finite joint support. Weights maximize
`ell(w) = sum_d log((A w)_d)` over the probability simplex, where `A` contains
paired donor likelihoods. Deterministic vertex exchange uses an exact line
search. Negative donor association is allowed; no covariance clipping or
arbitrary positive-variance floor is applied.

For `D` donors, `g = A' * (1 / (A w))`. Concavity gives
`(ell(w*) - ell(w)) / D <= max(g)/D - 1` on this fixed support, since `w' g = D`.
The reported tolerance is `1e-8`, with at most 20,000 iterations. Nonconvergence
is retained explicitly and cannot be used for prediction. This certificate
bounds the fixed-grid objective only; it does not prove continuous-support
optimality, a uniquely identified mixing distribution or biological calibration.
The general mixture-likelihood approach is described by
[Koenker and Gu, REBayes (2017)](https://www.jstatsoft.org/article/view/v082i08);
this implementation does not call REBayes.

Control-only query likelihood reweights the frozen joint support in log space.
The output includes latent control, treated and response moments, log1p-rate
moments, and the Poisson, cell-overdispersion and latent-rate components of
planned treated-count variance. Cell dispersion and fitted mixing parameters
are treated as fixed. Their estimation uncertainty is not integrated.

## Results on 2026-09-11

| Check | Kang | HIRISA |
| --- | ---: | ---: |
| Training cells | 2,651 | 119,513 |
| Exact depth groups | 2,375 | 50,251 |
| Available paired genes / requested | 8 / 16 | 11 / 16 |
| Converged finite-grid fits across four grids | 32 | 44 |
| Unavailable fits across four grids | 32 | 20 |
| Available control-only predictions across four grids | 1,984 | 2,728 |
| Unavailable predictions across four grids | 1,984 | 1,240 |
| Reconstructed training/query histograms | 1,248 | 1,152 |
| Independently compared numerical values | 3,347,684 | 4,461,776 |
| Maximum scaled numerical discrepancy | 1.33e-13 | 6.79e-13 |
| Grid refinement pass / available genes | 0 / 8 | 2 / 11 |
| Native peak RSS | 20,414,464 bytes | 111,689,728 bytes |
| Native wall time including stream handling | 1.57 seconds | 3.05 seconds |

At grid 65, 1,178 predictions are available and 806 unavailable. Four resolutions
of the same donor/gene combinations are not additional independent samples.
Kang's unavailable cases lack enough within-donor count pairs; HIRISA's five
unavailable paired cases have an unsupported control dispersion. No replacements
were invented. The exact gene statuses and all fits are in the retained reports.

`verify.py` reconstructs depth histograms and rate likelihoods from retained
individual-cell counts, independently solves a convex dual with one variable
per donor using SciPy SLSQP, and checks the objective, directional-gradient
certificate, all weights, latent moments and query/count-sampling arithmetic.
The dual upper bound is valid without assuming that learned mixing weights are
unique. All 76 available fits pass; there are no independent numerical errors.
Forty Swift regression tests in six suites pass on a physical Apple M4 Pro.
This scoped build does not qualify the complete Swift package or Metal execution.

The refinement protocol was frozen before fitting. Grid 33 to 65 must change
mean donor log likelihood by at most `1e-5` and every specified query moment by
at most `0.01 * max(1, abs(fine))`. Both conditions must pass. Only HIRISA MZB1
and PPP1R18 pass; 17 available models fail and 13 are unavailable. Kang's largest
scaled query-moment change is 0.432775; HIRISA's is 0.039998. HIRISA's largest
mean donor log-likelihood change is 0.028037. The tolerance was not relaxed after
seeing these failures.

Next: establish support convergence, account for fitted-parameter uncertainty,
qualify all genes and donor exclusions, then freeze and evaluate biological
outcomes on new data. The previous HIRISA-to-GSE181897 nominal 95% treated
coverage of **35.54% remains a failure**. These new conditional calculations
neither score nor repair it.

## Retained failures and numerical repairs

The first two test runs failed because Swift selected `Double.init(bitPattern:)`
for `UInt64.map(Double.init)`. Explicit `Double($0)` conversion repaired the
rate-likelihood score. Those source versions and failing logs are retained; no
real-data fit used the erroneous conversion.

The first complete real-data runs then exposed ten Kang query predictions whose
constant treated support had tiny positive variance around a rounded mean.
Constant-support axes now use their exact coordinate mean. Mathematical point
masses are classified before exponentiating log weights; underflowed components
and numerical variance collapse are reported separately. All original fits and
queries, the original executable and the comparison against corrected replay
are retained. Rates, fitted weights, likelihoods and all query probabilities are
exactly unchanged by this repair. The final model reports 220 Kang point-mass
predictions across the four grids and none for HIRISA. A fitted point mass is
not a claim that parameter or biological uncertainty is zero.

## Reproduction and retention

The [evidence manifest](evidence/2026-09-11/manifest.json) binds every logical file,
unique object and archive part. [Results](evidence/2026-09-11/results.json) retain
machine-readable totals and refinement failures. Restore into a new directory:

```sh
python3 Tools/Omics/CountObservation/Joint/retain.py restore \
  Tools/Omics/CountObservation/Joint/evidence/2026-09-11 /tmp/joint-evidence
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 \
  /tmp/joint-evidence/recipes/verify.py /tmp/joint-evidence/study Kang
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 \
  /tmp/joint-evidence/recipes/verify.py /tmp/joint-evidence/study HIRISA
```

Python requires NumPy and SciPy; preparation additionally needs h5py. To rebuild
on macOS with the declared scoped source dependencies:

```sh
bash Tools/Omics/H5AD/build.sh /tmp/joint-build
bash Tools/Omics/CountObservation/Joint/test.sh /tmp/joint-build
gzip -dc /tmp/joint-evidence/study/Kang-input.json.gz | /tmp/joint-build/joint-counts
```

`prepare.py STUDY Kang` and `prepare.py STUDY HIRISA` describe the original sparse
source replay, with the parent calibration and source paths declared in that
recipe. `run.py` is the pinned SSH execution recipe used in this study. Standalone
independent verification needs only the restored files, not the full H5ADs or
parent study directory. Retention includes individual-cell panel counts and full
depths, original-query cells, source hashes, all fits/queries, independent duals,
protocols, source snapshots, exact final/initial binaries and all failure logs.

The subsequent [adaptive support qualification](Adaptive/README.md) resolves the
tested support-convergence problem for all 19 available models from two initial
grids. Original grid-refinement failures above remain retained; parameter
uncertainty and biological validation are still separate open gates.
