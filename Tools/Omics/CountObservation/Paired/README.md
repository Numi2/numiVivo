# Paired count moments before treatment-response integration

The fitted [count observation marginals](../Calibration/README.md) cannot be
connected to treatment-response uncertainty by independently subtracting cell
measurement noise and using the remaining donor covariance. On all original
training genes, that calculation is indefinite for **7,489 Kang genes and 6,373
HIRISA genes**. A covariance matrix with a negative eigenvalue cannot describe
any joint distribution. Clipping the two marginal variances separately to zero
does not repair these failures.

`VivoPairedCountMoments.analyze` is a native, source-bound diagnostic for this
specific integration problem. It keeps complete paired donor counts and moments,
checks donor/condition identities, recomputes the published marginal calibration,
and reports the raw joint covariance, its minimum eigenvalue and the response
variance. It also compares the per-cell-rate and pseudobulk endpoints. It does
**not** fit a joint likelihood, project the covariance, produce treatment
predictions or replace the previous intervals.

## Estimand and calculation

For donor `d` and condition `c`, let `m_dc` be the average of the original
per-cell normalized rates `Y_i / (L_i / 1e6)`, including zero counts. The fitted
measurement variance is `v_dc = (P_dc + phi_c * R2_dc) / N_dc`, using all source
RNA features in each `L_i`. Every donor has equal weight. Control and treated
strata are paired by explicit donor identity, never array position.

Let `S` be the sample covariance matrix of paired `(m_d0, m_d1)` over donors.
Under independent measurement errors in the two distinct cell samples, the
method-of-moments latent covariance estimate is

```text
A = S - diag(mean_d(v_d0), mean_d(v_d1))
Var(latent treatment response) = A00 + A11 - 2*A01
```

The implementation calculates observed response variance directly from paired
differences, then subtracts both measurement variances. The covariance correction
does not assume that control and treatment donor rates are independent; it
assumes their cell measurement errors have zero covariance. Correlated technical
errors would require an additional observation model.

Both raw diagonal values are retained, including negative values. A second
diagnostic replaces only those diagonal values by their nonnegative parts,
matching the separate marginal-boundary operation in the published calibration.
Neither matrix is returned as a fitted probabilistic model. Minimum eigenvalues
use a scaled symmetric 2-by-2 formula; values below `-1e-10 * max(abs(entries))`
are classified as indefinite. This arithmetic tolerance is not variance
regularization. Zero matrices are positive semidefinite. The independent
reference uses NumPy/LAPACK's symmetric eigensolver.

If a condition's cell dispersion is unavailable, the joint correction remains
unavailable, with both marginal reasons retained. A Gamma-prior range failure
alone does not prevent diagnosing otherwise available measurement moments.
Negative donor covariance is retained; it is not itself an invalid covariance.

The other endpoint is RNA-weighted pseudobulk CPM, `1e6 * sum(Y_i)/sum(L_i)`.
It differs from the average of per-cell CPM when cell depth and expression rate
vary together. The existing pseudobulk response uses the donor-average difference
of `log1p(pseudobulk CPM)`. The new diagnostic compares that with the donor-average
difference of `log1p(mean per-cell CPM)`. Neither is the mean of individual-cell
log expression. The full feature vectors and maximum per-stratum differences
are retained, without filtering small or inconvenient effects.

## Complete original-data results

The input is the previously verified full set of sufficient moments from
**122,164 cells, 13 donors and 26 donor-condition groups**. Original H5AD cell
streams are not repeated. Their hashes, metadata and count bindings are inherited
explicitly from the completed calibration experiment.

| Full training fit | Kang | HIRISA |
| --- | ---: | ---: |
| Donors | 8 | 5 |
| Full RNA features | 15,706 | 18,082 |
| Both marginal calibrations available | 8,414 | 13,267 |
| Indefinite raw covariance | **7,489** | **6,373** |
| Positive semidefinite raw covariance | 925 | 6,894 |
| Joint correction unavailable | 7,292 | 4,815 |
| Negative noise-corrected response variance | 4,690 | 3,271 |
| Mean log-response direction differs between endpoints | 699 | 1,537 |
| Median absolute log-response endpoint difference | 0.04356 | 0.01561 |
| Largest absolute log-response endpoint difference | 0.52347 | 0.80999 |

Separate marginal clipping leaves exactly the same covariance classifications
on these full fits. Endpoint direction differences are descriptive arithmetic
sign changes, **not** statistical significance or evidence that either endpoint
is biologically preferable. Positive semidefiniteness is necessary for a joint
covariance, but does not prove positive-rate support or a calibrated model.

Every donor is omitted once and the calibration is recomputed over all remaining
donors and all original features: eight Kang and five HIRISA sensitivity fits.
These are not held-out predictions; no omitted donor outcome is scored.

| Across full fit and all donor omissions | Kang | HIRISA |
| --- | ---: | ---: |
| Genes whose covariance changes sign | 4,786 | 3,392 |
| Genes whose covariance status changes | 3,093 | 8,264 |

Status changes include becoming available or unavailable, as well as switching
between indefinite and positive semidefinite. Counts cover all original genes,
including those with unavailable marginal calibration. Every omitted-donor
report and every feature's sensitivity record are retained.

These findings are compatible with limited donor information, unstable
measurement-noise estimation and violated model assumptions. They do not
identify a unique biological cause. They specifically rule out treating the
separately corrected moments as an automatically valid joint donor model.

## Execution and numerical verification

The physical Apple M4 Pro runs the actual NumiVivoKit implementation, built with
Swift 6.3.3. **34 tests in five suites pass**, including explicit donor pairing,
missing/duplicate-stratum rejection, a known valid covariance, an impossible
covariance despite positive diagonals, negative association, full-library
endpoint arithmetic and absent count support.

All **15 reports and 249,846 feature records** are independently reconstructed
from the original NumPy cell moments. The reference checks 4,106,588 finite
numerical values, every feature and donor identity, source fingerprints, all
missing-value masks and all calibration/covariance statuses. Maximum error is
**4.35e-11** on `abs(native-reference)/max(1,abs(reference))`, below 2e-7.

The complete native output is 138,878,820 bytes for Kang and 120,786,013 bytes for
HIRISA before lossless compression. Peak native RSS is 192,643,072 and
172,949,504 bytes respectively. Streamed times are 5.37 and 8.12 seconds; these
include SSH/output transport and are not controlled CPU/Metal speed comparisons.
The exact executable is hashed before both runs and retained. Its bytes remain
identical after the final test-harness relink, which happened only after both
scientific runs finished. No fresh full-source replay was needed for this
moment-only analysis.

An initial source-transfer attempt failed because the remote Python version
lacked a requested tar extraction argument; the build then found no source.
The corrected transfer validated regular relative members before extraction
into the still-empty directory. Both admission failures are retained. Existing
unrelated compiler warnings remain in the build log. The initial 33-test run
and final 34-test run are both retained.

## Reproduction and evidence

```sh
bash Tools/Omics/H5AD/build.sh /new/build
bash Tools/Omics/CountObservation/Paired/test.sh /new/build
gzip -dc /restored/study/Kang-input.json.gz | /new/build/paired-counts
```

The executable consumes the frozen JSON on standard input and emits one JSON
line per full or donor-omission result. `prepare.py` freezes the complete parent
moments; `run.py` is the physical-Mac SSH adapter with a pinned binary identity.
The native library interface remains the owner; these are qualification tools,
not newly advertised production CLI commands.

The [manifest](evidence/2026-09-11/manifest.json) retains full native outputs,
full independent references, source moments, input identities, all sensitivity
vectors, exact executable, source snapshots, logs and failures. The gzip tar is
split into three lossless parts below repository file-size limits; the manifest
checks each part and the joined archive before extraction. All 51 logical files
were restored and checked. Restore and
verify both origins without fetching the original H5ADs:

```sh
python Tools/Omics/CountObservation/Paired/retain.py restore \
  Tools/Omics/CountObservation/Paired/evidence/2026-09-11 /restored
python /restored/recipes/verify.py /restored/study Kang /restored/reference
python /restored/recipes/verify.py /restored/study HIRISA /restored/reference
```

Original H5ADs and raw cell-stream admission remain bound to the separately
published [calibration evidence](../Calibration/evidence/2026-09-11/manifest.json).
This archive is sufficient for the paired-moment computation, not a substitute
for the original source-data qualification.

## Next model requirement

Fit a coherent joint latent donor-response model with explicit positive-rate
support and cell measurement likelihoods, including negative donor association
where the data support it. Estimate covariance jointly with observation noise;
assess dispersion shrinkage and parameter uncertainty rather than inserting
the invalid matrix or declaring a PSD projection biologically calibrated.
Declare the per-cell or pseudobulk endpoint before scoring. Integrate only query
control information and explicit planned sampling for future treated counts.

The earlier GSE181897 point-prediction results and failed treatment interval
coverage remain unchanged. Known training/development data and a valid numerical
implementation do not replace new independent biological validation. The wider
twelve-part development goal remains active.

The subsequent [joint count-response model](../Joint/README.md) now fits paired
nonnegative rate distributions on the unchanged 16-gene development panel.
Finite-grid numerical checks pass, but 17/19 available models fail grid
refinement; this does not repair or supersede the full-gene diagnosis above.
