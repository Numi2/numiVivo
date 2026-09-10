# Native global QL scale and refit

`VivoOmicsNBQLGlobalScale.fit` now estimates the global quasi-likelihood scale
for an explicit pseudobulk feature family. It uses the existing native NB
solver, the [adaptive conditional moment owner](../Moments/Adaptive/README.md)
and a new native robust LOWESS implementation. Counts, donor design, offsets,
contrast, positive trend dispersions and abundance coordinates are supplied by
the caller. This is a CPU FP64 stage; existing cohort Wald/LRT defaults remain
unchanged.

The method fits counts at the supplied trend dispersions, performs exactly two
scale updates at those fixed initial means, then refits once at trend dispersion
divided by the final scale. Each update smooths the fourth root of adjusted
quasi-dispersion over abundance and takes the type-7 90th percentile of the
smooth, floored at one and raised to the fourth power. Final adjusted residuals
use the refitted means. Abundance coordinates in this qualification come from
the pinned reference; native abundance estimation is not qualified here.
The subsequent [native abundance integration](../Abundance/README.md) qualifies
that missing input through an explicit native entry point on all 58 arms.

`VivoOmicsRobustLowess.fit` uses tricube neighborhoods, stable local linear
regression, three Tukey-bisquare robustness updates and delta interpolation.
Default span is 0.5 and delta is 0.01 times the abundance range. Sorting preserves
input order for ties, and results return in input order. Constant fits,
empty-weight anchors, robustness iterations and neighborhood work are explicit
diagnostics. Global scale estimation refuses empty weighted neighborhoods.
The neighborhood counter reserves two passes per visited interval; it is a
work-budget counter, not a throughput metric. The method follows standard
[LOWESS](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/lowess.html)
and the global-scale stage of [edgeR v4](https://doi.org/10.1093/nar/gkaf018),
without importing reference package source into the native owner.

Initial or refit errors and unidentified/unconverged genes produce a structured
incomplete family with retained results and feature/stage failures. Unavailable
quasi-dispersions are explicitly omitted from scale learning; fewer than three
informative residuals fails the update. No gene is silently dropped to pass a
family. Cancellation propagates. The public bounds are 3–100,000 genes and at
most one million gene-by-donor count entries, subject to the existing QR donor
bound. Defaults allow 100 fit iterations, one million moment evaluations per
observation and 100 million neighborhood visits per smoother. Unrepresentable
scaled dispersion or exhausted work returns a failure, not substituted evidence.

## Qualification

The [frozen protocol](PROTOCOL.md) retains all 29 Kang, Hagai and Crowell cases
and both original full-support trend arms. All 58 native arms complete. All 58
independent checks pass across 483,576 gene/arm fits, each with an initial fit
and a final refit. Neither scale update omits any gene in this family.

| Check | Observed maximum | Frozen limit |
| --- | ---: | ---: |
| Controlled native LOWESS vs R, 25 cases | 7.02e-10 | 2e-7 |
| LOWESS on all 116 actual native update inputs | 1.05e-8 | 2e-7 |
| Initial/final means and contrast effects vs independent fixed-scale GLMs | 1.93e-6 | 2e-5 |
| Native scaled score, independently recomputed from counts and means | 9.99999e-8 | 1e-7 |

Relative errors in this table divide by `max(1, abs(reference))`. Independent
score evaluation uses the NB score and Fisher diagonal directly; the checker
allows 1e-10 absolute rounding discrepancy, and the observed discrepancy from
the native reported score is zero. All 31 focused Swift tests in six suites
pass on the physical M4 Pro with Apple Swift 6.3.3. Tests cover the scale floor,
fixed initial means between updates, scaled refits, family failures and work
limits as well as existing NB, cohort, support and moment behavior.

The independent GLMs call `edgeR::mglmLevenberg` directly with tolerance 1e-12
and at most 10,000 iterations, using the exact native final scale. Native-input
LOWESS checks use `stats::lowess`. Reference versions are edgeR 4.10.5,
limma 3.68.5 and jsonlite 2.0.0. Reference extraction reconstructs all original
29-case scales/refits and checks their hashes before use.

The largest relative global-scale difference from the original edgeR fits is
2.33e-5. Maximum final adjusted residual-deviance and residual-DF differences,
divided by `max(1, abs(reference))`, are 0.00122 and 0.000610. These are
descriptive comparisons: the moment approximation and original reference fit
stopping behavior differ. They are not failures of an asserted edgeR identity
or evidence of better statistical calibration. Full differences remain in the
per-gene checks.

The run used 5,318,700,013 direct-count or adaptive-endpoint evaluations across
the three residual stages. Observed arm wall times range from 5.67 to 76.78
seconds, including serialization, SSH and shared-host effects. These are not a
controlled speed comparison. The native executable SHA256 is
`d9d5a8b989cb9968792caf1185ae9911634f8e4872865af815495282b93b3fd6`.

## Retained failures and remaining work

The first test compilation could not type-check a mixed UInt64/integer literal
expression. Explicit types repaired the test fixture; production source and
the qualified binary did not change. The first independent R checker collapsed
JSON updates into a dataframe. Disabling that simplification repaired parsing.
The second checker passed `tol`/`maxit` through `glmFit`, whose pinned wrapper
ignores those arguments. Its comparisons failed at default solver precision.
Calling `mglmLevenberg` directly applies the intended limits and all final
comparisons pass. Both failed checker attempts, their outputs/logs and the
initial test compilation remain archived.

Robust abundance-dependent prior estimation with unequal residual DF,
posterior moderation, constrained QL hypothesis testing and cohort integration
remain open. This stage does not supply those methods or establish FDR,
coverage, power, biological truth, varying-support borrowing, million-cell
execution or GPU acceleration. In particular, the
parent study's 39 Hagai sham calls with adjusted native-trend reference QL
versus 20 with Wald remain a calibration concern.

## Reproduce and inspect

```sh
bash build.sh NATIVE_RUN/build
python check_lowess.py --out RUN/lowess --binary NATIVE_RUN/lowess
python reference.py --ql-root PRIOR_QL_RUN --out RUN/reference --jobs 1
python run.py --ql-root PRIOR_QL_RUN --reference-root RUN/reference --out RUN/native --binary NATIVE_RUN/build/nb-ql-global --jobs 3
python check_reference.py --ql-root PRIOR_QL_RUN --reference-root RUN/reference --root RUN/native --jobs 1
python summarize.py --root RUN --ql-root PRIOR_QL_RUN
python archive.py --root RUN --out ARCHIVE
python archive.py --out ARCHIVE --verify
```

The standalone LOWESS probe is compiled with Swift 6 from
`VivoOmicsLinearStatistics.swift`, `VivoOmicsRobustLowess.swift` and `Lowess.swift`
using `swiftc -O -parse-as-library`; `build.sh` compiles the full global stage.
Reference scripts use the pinned R library and explicit prior study paths.
The archive step also requires `execution.json` and collected native build/test
logs; its source hashes must match the tested checkout. Existing run receipts
are reused only after payload, executable and output hash checks.

The [summary](evidence/2026-09-10/summary.json.gz) and
[manifest](evidence/2026-09-10/manifest.json) bind complete receipts, per-gene
fit/score/residual comparisons, native update diagnostics, smoother references,
controlled cases, source/test hashes and retained failures. Large native,
reference and independent matrices remain external at their exact listed
hashes. The archive verifier checks stored members; it does not claim that
external files are available on an arbitrary reader's machine.
