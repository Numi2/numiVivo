# Adaptive native deviance moments

`VivoOmicsNBResidualAdjustment` accepts an explicit `method: .adaptive` for
`moments` and `adjustedResiduals`. The existing default remains `.direct`, and
its published work-limit failures remain in the [parent evidence](../README.md).
This method changes how the conditional NB moments are summed. It does not
estimate global QL scale, robust priors or cohort hypothesis tests.

For NB means at least 10,000, adaptive blocks sum a cubic Hermite interpolant
of probability, probability times deviance, and probability times deviance
squared. A fourth-derivative remainder bounds each finite integer block.
Smaller means and exact Poisson inputs retain direct summation. The
[derivation](DERIVATION.md) gives the polynomial identity, derivative bounds and
normalization error propagation. The [protocol](PROTOCOL.md) was frozen before
native adaptive results.

Success requires both interpolation and omitted-tail errors to satisfy the
requested tolerance. `summation` records the method, covered integer support,
accepted blocks and raw/normalized interpolation bounds. `evaluatedCounts`
counts actual endpoint evaluations on the adaptive path. The mean/variance
summation and truncation bounds must be added; the combined relative bound is
at most 1e-10 under the default settings. Floating-point rounding is separate
from these analytic bounds. Work exhaustion returns an error.

The original twenty resource-limited gene/arm fits now all complete within the
unchanged one-million-evaluation budget. Independent bounded-chunk PMF summation
checks all 124 moments, including the six previously unavailable observations.
These fits use 4,455,408 endpoint/direct evaluations in total, with at most
43,482 for any observation. Comparison with the 118 previously available
higher-work direct moments is retained. All twelve extended-grid cases pass;
maximum relative error across those and the recovered moments is 2.82e-9.

All 58 arms across the original 29 Kang, Hagai and Crowell cases pass, covering
3,940,972 moments. Of these, 24,288 use adaptive blocks and the remainder retain
direct summation. Every previously available moment is compared with the
retained direct result; maximum relative mean/variance/scale/DF difference is
1.10e-11. The newly available 124 moments match the separately checked twenty-fit
outputs exactly. Counts, fitted means, trend dispersions, supplied scales,
designs and feature families retain their original hashes. All combined
relative approximation bounds are at most 1e-10. This stage used 1,769,552,309
direct-count or adaptive-endpoint evaluations; those operations have different
costs, so this is a work inventory rather than a throughput comparison.

All 28 focused Swift tests in five suites pass on the physical M4 Pro, including
polynomial identities, exact quartic remainders, independent 65-digit PMF/slope
and fourth-derivative checks, direct/adaptive comparisons and work-budget errors.
The original 77-case grid also passes. Native computation is CPU FP64. The
executable SHA256 is
`c4759843fbe6e061621b2d2d032892c04ea114ab3e94a71c1ef85e9bccb1820e`.

The first pilot's native-trend arm took 170.64 seconds while the second arm took
4.18 seconds. A retained native process sample found the first waiting inside
output writing. Shared-host transport and JSON overhead are included in these
measurements. This is not a controlled end-to-end speed comparison or a claim
about million-cell execution, GPU acceleration, calibration or biological truth.

## Reproduce

```sh
bash build.sh NATIVE_RUN/build
python check.py --root RUN --old PRIOR_DIRECT_RUN --binary NATIVE_RUN/build/nb-adaptive-moments
python check_grid.py --out RUN/grid --binary NATIVE_RUN/build/nb-adaptive-moments
python run.py --ql-root PRIOR_QL_RUN --out RUN/family --binary NATIVE_RUN/build/nb-adaptive-moments --jobs 3
python compare_family.py --root RUN --old PRIOR_DIRECT_RUN --ql-root PRIOR_QL_RUN
```

`check.py` expects `twenty-native.json` generated from `twenty-input.json` by the
native executable. That input wraps the exact prior work-sensitivity rows in a
`residuals` object. Existing completed outputs are not overwritten by a rerun.
Complete-family coverage requires all 58 arms; a pilot alone does not qualify it.

The [archive manifest](evidence/2026-09-10/manifest.json) binds complete receipts,
per-gene comparisons, grids, recovered results, native source/test hashes and
logs. Large family arrays remain external at exact listed hashes. Initial
compiler warnings, a failed expression-format repair and the test-fixture
indentation failure are retained; they precede the final successful checks.
The [family comparison](evidence/2026-09-10/family-comparison.json.gz) records the
complete coverage. The direct method's original failures remain in the parent
archive. The subsequent [native global scale/refit](../../GlobalScale/README.md)
passes the complete 58-arm family using explicit abundance inputs. Robust
unequal-DF moderation and cohort inference remain open, as do native abundance
estimation and biological/FDR calibration.
