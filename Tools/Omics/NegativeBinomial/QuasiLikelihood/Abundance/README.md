# Native NB abundance and global QL integration

`VivoOmicsNBAbundance.fit` estimates average log2 counts per million from exact
counts, log library offsets and a supplied scalar NB dispersion. The explicit
`VivoOmicsNBQLGlobalScale.fitEstimatingAbundance` entry point uses these native
covariates for the [global QL scale/refit stage](../GlobalScale/README.md).
Reference abundance values are no longer required by this path. The existing
supplied-covariate entry point and cohort Wald/LRT defaults remain available.

For library size `L_i = exp(offset_i)` and average prior count `p` (default 2),
the abundance response is `z_i = count_i + p*L_i/mean(L)`. Its offset library is
`L'_i = L_i + 2*p*L_i/mean(L)`. An intercept-only NB2 fit has means
`mu_i = L'_i*exp(beta)` and score

```text
S(beta) = sum_i (z_i - mu_i) / (1 + phi*mu_i)
S'(beta) = -sum_i mu_i*(1 + phi*z_i) / (1 + phi*mu_i)^2
log2CPM = (beta + log(1,000,000)) / log(2)
```

The monotone score has a unique finite root when the augmented total is positive.
This follows the abundance definition documented by
[edgeR](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf)
(`aveLogCPM`, `addPriorCount` and `mglmOneGroup`). The native implementation uses
the mathematical score rather than reference package code.

Offsets are centered before fitting. Prior scaling uses centered exponentials;
the common offset adjustment is evaluated in log space. The positive-dispersion score uses
bounded logistic factors. A sign bracket followed by bisection stops at a
computed bracket width of at most 2e-12, or an exactly zero floating-point score.
Floating-point rounding is distinct from this computed bracket width. The
returned result retains log proportion, log2 CPM, dispersion, prior count,
score evaluations, bracket width and scaled score. Phi zero uses the exact
Poisson score. No continuous pseudo-count is passed to the original count GLMs.

The scalar interface accepts 1–512 observations, raw counts at most 2^53,
finite offsets in [-700,700], phi zero or [1e-8,100], and prior counts in
[0,1e8]. The default budget is 128 score evaluations; callers may choose
1–1,024. A positive scaled prior that underflows is rejected. With prior zero,
an all-zero row has no finite log abundance and returns an error. Positive-prior
zero rows have finite abundance, but still fail the existing count GLM if their
positive-count support is unidentified. Observation weights and gene-by-donor
offset matrices are outside this scalar/common-offset interface.

The integrated entry point preserves each abundance fit and indexed failures.
An abundance failure leaves the global fit unavailable; a later count-fit or
scale failure retains both the abundance results and the global stage's partial
results. Cancellation propagates. Original counts, designs, contrasts and trend
dispersions enter the NB family unchanged. The family retains the existing
3–100,000 gene, at-most-one-million gene-by-donor-entry bounds.

## Controlled qualification

All 105 controlled cases pass: 100 finite comparisons and five expected errors.
The cases cover zero counts, prior counts 0/0.5/2/10, Poisson and NB dispersions,
equal and unequal libraries, single observations, exact counts up to 2^53,
wide/small/large libraries, count representability and exhausted work. Every
finite result agrees with pinned edgeR 4.10.5 and an independent 80-digit score
root. Maximum native-versus-high-precision absolute error is 1.32e-12 log2 CPM,
below the frozen 2e-8 limit. There are no reference disagreements.

All 35 focused Swift tests in seven suites pass on the physical M4 Pro with
Apple Swift 6.3.3. New tests also verify offset-shift and permutation identities,
closed-form equal-library fits, the unchanged original count response, and
retention of abundance versus count-fit failures. Native computation is CPU
FP64. The executable SHA256 is
`05a2875637aa4626c5a5ee39a7b567aff913be4905ad399df12d1dfccbcc7a76`.

## Complete real-data qualification

Both original trend arms pass for all 29 Kang, Hagai and Crowell cases: 58 arms
and 483,576 gene/arm fits, with no gene omitted from either global-scale update.
Every initial NB fit is identical to the previously qualified supplied-abundance
result. All abundance score roots are independently bracketed from reconstructed
augmented counts/libraries, and all final count-model scores are independently
recomputed from original counts and the new fitted means.

| Comparison | Maximum observed error | Frozen limit |
| --- | ---: | ---: |
| Native abundance vs pinned reference covariate, absolute log2 CPM | 2.20e-11 | 2e-8 |
| Global scale vs prior native supplied-abundance fit, relative | 6.18e-13 | 2e-7 |
| Refit means vs prior native supplied-abundance fit, relative | 2.39e-13 | 2e-7 |
| Final adjusted deviance/DF/quasi-dispersion vs prior native fit, relative | 5.82e-13 | 2e-7 |

Relative errors divide by `max(1, abs(reference))`. Maximum independently
computed abundance scaled score is 1.81e-10; the final count-model maximum is
9.99999e-8, within the native 1e-7 limit. Independent score checks allow 1e-10
absolute rounding discrepancy. Abundance estimation uses 20,310,187 score
evaluations across this family, at most 42 per gene/arm. This is a work inventory,
not a controlled speed comparison.

One original transfer stalled after the native result reached output writing.
A native stack sample captured `write` beneath `NSConcreteFileHandle.writeData`,
and socket observations retained the queued output. After 356.50 seconds, only
that SSH transfer was terminated; its exit-255 receipt, partial output and
original incomplete family remain preserved. The native-trend Crowell-07 arm
was rerun with the unchanged binary, saving output on the Mac mini before a
compressed IPv4 transfer. It completed and passed in 11.70 seconds including
transport. The other 57 results were reused unchanged. This demonstrates the
recovery path, not general throughput superiority.
All 39,780,352 bytes captured before the failed transfer match the prefix of
the recovered 63,021,046-byte native result; the prefix check is retained.

An earlier process sample missed its already-completed target. Local sampling
attempts also failed: the PATH command resolved to a broken Python program, and
explicit system samplers did not return and were stopped. No local stack claim
is made from those attempts; their logs are retained separately from the
successful native stack sample.

The [summary](evidence/2026-09-10/summary.json.gz) and
[manifest](evidence/2026-09-10/manifest.json) bind complete per-gene comparisons,
controls, receipts, source/test hashes, native logs and the retained transport
failure. Large full and partial native outputs remain external at exact hashes;
the manifest distinguishes them. The archive verifier checks stored members,
not the availability of external files on another machine.

## Reproduce

```sh
bash build.sh NATIVE_RUN/build
python check.py --out RUN/controls --binary NATIVE_RUN/build/nb-ql-abundance
python run.py --ql-root PRIOR_QL_RUN --prior-global-root PRIOR_GLOBAL_RUN --out RUN/family --binary NATIVE_RUN/build/nb-ql-abundance --case kang-null-01-default
python run.py --ql-root PRIOR_QL_RUN --prior-global-root PRIOR_GLOBAL_RUN --out RUN/family --binary NATIVE_RUN/build/nb-ql-abundance --remote-spool-root NATIVE_RUN/spool --jobs 3
python summarize.py --root RUN
python archive.py --root RUN --out ARCHIVE
python archive.py --out ARCHIVE --verify
```

The reference check pins edgeR 4.10.5, limma 3.68.5 and jsonlite 2.0.0, and
uses mpmath at 80 decimal digits for independent roots. The native harness runs
on `macmini`; Python and R are qualification tools. Before archiving, collect
native build/test logs and `execution.json` with source, executable and prior
artifact hashes. Completed family receipts are reused only after input,
executable, prior-output and result hash checks. Controlled reruns require a
new output directory so prior attempts remain intact.
`--remote-spool-root` writes native results without overwriting existing remote
files and retains a local execution receipt before fetching compressed output.
A failed download can then retry the transfer without rerunning computation.
Keep that remote directory and execution receipt until the output is verified.

The [protocol](PROTOCOL.md) fixes all 58 original full-support arms and their
acceptance limits. The subsequent [native robust moderation](../Moderation/README.md)
implements unequal-DF prior fitting and posterior variance, with all 58 arms
qualified. Cohort QL hypothesis testing remains open. This stage does not establish FDR,
coverage, power, biological truth, varying-support borrowing, general speed
superiority, GPU acceleration or million-cell execution. The parent study's
39 Hagai sham calls with adjusted native-trend reference QL versus 20 with Wald
remain a calibration concern.
