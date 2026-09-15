# Kang B-cell Metal NB execution attribution (2026-09-15)

This record adds a bounded execution receipt for the experimental Metal FP32
negative-binomial objective. It does not promote that path or change the
CPU-authoritative policy recorded in the parent Kang evidence.

## Bound invocation

- Implementation source: `f5fab03c65057705a401e9cd3ca8bb6c26bd0684`
  (`Profile opt-in Metal NB execution`).
- Prepared Kang B-cell H5AD:
  `abbc78e383b298ba8634dba0f267219f75092dfaff54097bbe8e9f7d73391802`
  (2,651 cells, 15,706 features, 1,479,543 nonzeros, eight donors and 16
  donor-condition pseudobulks).
- Metal plan:
  `ef325faa1e0acce478d574f1e498f7e25291a0b2727d8129643f0d7606da8935`;
  it selected `negativeBinomial`, `metalFP32`, a parametric trend, and the
  paired-donor IFNB-versus-control contrast.
- The emitted sidecar and analysis receipt respectively hash to
  `c513877d63e3902d8431d31908eceefa0a0c237499f1086c570b9170211bd7f3`
  and `c1765dc5552dcc2f642a4525896bdc0cdb9892042d3f02afa0574d32c21e8dab`.

The input, native store and result payloads are local retained artifacts, not
committed data. The command completed the analysis in 225.02 seconds and
emitted the optional profile sidecar only after writing the ordinary analysis
receipt.

## Measured execution

The profile records a 221.561-second cohort evaluation. It attempted and
completed 240,785 fixed-dispersion fits and 1,921,588 Metal objectives. All
objectives used a 16-observation pseudobulk vector: 30,745,408 observations in
1,921,588 batches. The 196.935 seconds attributed to Metal objective calls
included 186.747 seconds in command completion (94.8%), 7.543 seconds encoding
and submitting commands, 0.565 seconds shared-buffer writes, 0.501 seconds CPU
term reduction, and 0.172 seconds FP64-to-FP32 conversion. The retained CPU
exact-likelihood calls contributed 0.435 seconds across 1,140,012 evaluations.

This is host-local wall-clock attribution, not a GPU-profiler trace or proof of
causality. CPU dispersion profiling, import/export, aggregation, and other
analysis work are outside those component totals. In particular, this record
does not claim useful end-to-end Metal performance.

## Coverage and replay boundary

Of 8,894 eligible features, the run retained 3,803 tested features, 3,494
rank-deficient-support rejections, and 1,597 numerical failures. The latter
all reported an unconverged NB dispersion-profile coefficient fit. No acceptance
criterion was loosened and no failed Metal fit was silently replaced by CPU.
The existing CPU Kang evidence records 5,400 tested features on the same
eligibility scope, so this invocation does not demonstrate numerical-coverage
recovery or promotion.

`singlecell-analysis-verify` replayed the persisted analytical receipt in
219.22 seconds and returned
`verified-native-analysis-not-biological-calibration` with 3,803 tested
features. Timing is intentionally excluded from the persisted analytical result
because wall-clock values cannot survive source reconstruction; the sidecar
binds it to that receipt instead. Replay therefore validates analytical-content
reconstruction, not timing, performance, numerical agreement, statistical
calibration, or biological inference.

The result is a concrete reason to keep the Metal path experimental and to
investigate CPU optimization or batched resident-data execution before expanding
it. A current CPU timing comparison and CPU-versus-Metal agreement/coverage gate
remain required before any promotion.
