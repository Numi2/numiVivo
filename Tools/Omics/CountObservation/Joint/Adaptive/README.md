# Adaptive joint count-response support

This extends the [finite-grid model](../README.md) by searching the continuous
nonnegative control/treated rate domain for likelihood-improving support points.
It preserves fixed NB2 cell dispersions, the original training/query separation,
and the latent-rate endpoint. A small likelihood gap does not establish unique
mixing weights, stable posterior moments or calibrated biological uncertainty.

## Native method and domain bound

`VivoAdaptiveJointCountResponse` initializes the declared grid with nine points
per axis plus every training-donor MLE. It adds the oracle's control and treated
coordinates to the Cartesian support, reuses the previous fitted weights, and
runs the shared finite-support vertex-exchange optimizer. Adaptive fits choose
the active weight exchange with the largest exact-line-search likelihood gain
toward the maximum-gradient component; the legacy fixed-grid API retains its
previous selection. That optimizer keeps
its `1e-8` mean likelihood-gap tolerance. The continuous target is `1e-6`, with
at most 64 support additions, 4,096 oracle leaves and 200,000 inner weight
iterations. Exhaustion is explicit and refuses adaptive prediction.

For fixed fitted donor likelihoods `s_d`, the continuous directional score is

```
G(rC,rT) / D = sum_d L_dC(rC) * L_dT(rT) / (D * s_d).
```

The maximum score minus one bounds the mean log-likelihood improvement from any
mixing distribution over the rate domain. It is sufficient to search the product
of training-donor MLE ranges: outside either range, moving that coordinate
inward increases every donor's conditional count likelihood.

For each condition choose `S = 1 / (1 + max_i(phi * exposure_i))` and write
`r = S * expm1(x)`. The rate-dependent NB2 log likelihood becomes

```
Y * log(expm1(x))
  - sum_i (y_i + 1/phi) * log(1 - a_i + a_i * exp(x)),
a_i = phi * exposure_i * S < 1.
```

This is concave in `x >= 0`; the Poisson limit is also concave. A tangent to each
paired log likelihood therefore upper-bounds it over a rectangle. The sum of
exponentials of those affine tangents is convex, so its maximum on a rectangle
occurs at a corner. Each box uses the smaller of this bound and the sum of
individual donor maxima obtained by clamping each donor's MLE to that box.
The native search splits the box with the largest upper bound. An evaluated
score above `1 + 0.9 * tolerance` supplies a new support point, reserving
headroom for the floating-point allowance; otherwise subdivision
continues until the global bound passes or a budget is exhausted.

Every retained leaf includes its binary subdivision path, coordinate bounds
and directional upper bound. Leaves cover the original compact rectangle.
The explicit FP64 allowance is `1e-9 * (1 + rawBound)`. This is an analytic
bound evaluated and independently checked in floating point, **not a formal
directed-rounding interval certificate**.

This is a native implementation of a mixture-likelihood optimization approach;
[Koenker and Gu (2017)](https://www.jstatsoft.org/article/view/v082i08) provide
background on nonparametric mixture estimation. This code does not call REBayes.

## Experimental scope

Both original sources contribute all 122,164 admitted training cells and their
full RNA denominators. The experiment uses the unchanged identifier-selected
16-gene panel, including all unavailable cases. These are eight Kang and five
HIRISA training donors, not 122,164 independent donor observations. Query inputs
are the original controls from 62 GSE181897 donors, with the same planned future
treated-cell depths. No query-treated outcomes enter fitting or scoring.

`verify.py` reconstructs likelihoods and coordinate derivatives from retained
individual-cell vectors. It independently checks every final leaf bound, proves
complete nonoverlapping domain coverage through the prefix-free binary paths,
and compares fitted likelihoods, support witnesses and all admitted query
probabilities and moments. `compare.py` measures differences from the previous
65-point grid; those differences are model diagnostics, not outcome errors.

The first budget used 20,000 inner weight iterations. All 18 other available
models reached the continuous bound, but Kang PNPLA4 stopped at a finite-support
gap of `1.5862807574151816e-6`. Its predictions were explicitly unavailable.
Complete initial fits, queries, independent checks, source and executable are
retained. A second run with 200,000 iterations still stopped PNPLA4 at
`6.393003479931991e-7`. Both exhausted runs are retained. A focused independent
diagnostic reached the `1e-8` target in 13 largest-gain exchanges from that state,
which motivated the native selection repair. The final iteration budget remains
200,000 and both likelihood tolerances are unchanged.

The final measured results and reproducible evidence are described below.

## Final results, 2026-09-11

All **19 available origin/gene models** reach the continuous mean likelihood-gap
bound of `1e-6`, from both grid-9 and grid-17 initializations. Thirteen of the
32 requested origin/gene cases remain unavailable because of the original
cell-dispersion admission rules; none were replaced or silently dropped.

| Numerical check | Kang | HIRISA |
| --- | ---: | ---: |
| Available models per initialization | 8 | 11 |
| Available control-only predictions per initialization | 496 | 682 |
| Unavailable predictions per initialization | 496 | 310 |
| Independently checked values, grid-9 initialization | 272,129 | 414,158 |
| Independently checked values, grid-17 initialization | 451,981 | 681,760 |
| Final leaf bounds checked, both initializations | 9,292 | 6,180 |
| Maximum scaled reference discrepancy | 7.82e-14 | 5.19e-12 |
| Maximum scaled query-moment difference between initializations | 1.65e-6 | 6.71e-9 |

All **1,820,028** comparisons and **15,472** final leaf checks pass. All 19
available models pass the predeclared 1% initialization-sensitivity criterion
for treated mean/variance, log1p response mean/variance and planned total-count
variance. The criterion and exact denominator were frozen before the alternate
initialization. This is evidence for these two starts at these tolerances, not
proof of unique mixing distributions or statistical parameter uncertainty.

Forty-three Swift regression tests pass on the physical M4 Pro. The shared
solver refactor also reproduces all **24,985,744 expanded bytes** of the previous
four-grid Kang/HIRISA native outputs exactly. The original fixed-grid failures
remain historical evidence; the adaptive method supplies a separate numerical
resolution of that support problem.

This remains a 16-gene development-panel experiment using all admitted training
cells. The [full-gene extension](Full/README.md) now completes Kang fitting and
independent numerical checks, with nine eligible solver-limit cases; HIRISA is
still running. Donor-exclusion sensitivity, mixing-distribution and dispersion
estimation uncertainty, and independent treated-outcome validation remain open. The historical nominal 95% treated coverage of **35.54%** is unchanged;
no treated outcomes were scored or used to recalibrate intervals here.

## Reproduction

The [manifest](evidence/2026-09-11/manifest.json) retains source snapshots,
all three primary attempts, initializations, original individual-cell panel
counts and full RNA depths, query cells, exact binaries, partition certificates,
references, protocols, logs and the unchanged fixed-grid replay. The
[machine-readable results](evidence/2026-09-11/results.json) preserve every status.

```sh
python3 Tools/Omics/CountObservation/Joint/Adaptive/retain.py restore \
  Tools/Omics/CountObservation/Joint/Adaptive/evidence/2026-09-11 /tmp/adaptive-evidence
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 \
  /tmp/adaptive-evidence/recipes/verify.py /tmp/adaptive-evidence/study Kang
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 \
  /tmp/adaptive-evidence/recipes/verify.py /tmp/adaptive-evidence/study HIRISA
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 \
  /tmp/adaptive-evidence/recipes/verify.py \
  /tmp/adaptive-evidence/study/initialization17 Kang /tmp/adaptive-evidence/study/parent
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 \
  /tmp/adaptive-evidence/recipes/verify.py \
  /tmp/adaptive-evidence/study/initialization17 HIRISA /tmp/adaptive-evidence/study/parent
python3 /tmp/adaptive-evidence/recipes/compare_initializations.py /tmp/adaptive-evidence/study
```

The reference requires NumPy and SciPy. Build on macOS with
`bash Tools/Omics/H5AD/build.sh /tmp/adaptive-build`, then
`bash Tools/Omics/CountObservation/Joint/Adaptive/test.sh /tmp/adaptive-build`.
Pipe a restored `ORIGIN-input.json.gz` through `gzip -dc` into
`/tmp/adaptive-build/adaptive-counts`. The alternate initialization uses the
retained `initialization17/Run.swift`, whose only harness change is the initial
grid value. This scoped build is not full-package, GPU or biological qualification.
