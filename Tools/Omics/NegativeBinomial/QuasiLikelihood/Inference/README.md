# Native adjusted QL cohort inference

The experimental `quasiLikelihoodAdjusted` option now runs the complete native
abundance, global scale/refit, adjusted residual, robust moderation and constrained
hypothesis chain through the existing pseudobulk expression owner. It uses the
cohort dispersion trend, not the gene-wise MAP dispersion. Existing Wald/LRT
behavior and defaults remain unchanged; the sham results below do not support
promoting QL as a calibration fix.

## Use and interpretation

For an existing `VivoOmicsExpressionContrast` request:

```swift
request.model = .negativeBinomial
var options = VivoOmicsNBCohortOptions()
options.testMethod = .quasiLikelihoodAdjusted
request.negativeBinomialOptions = options
let result = try VivoPseudobulkDifferentialExpression.run(dataset, contrast: request)
```

The same option is available in the existing serialized analysis request as
`negativeBinomialOptions.testMethod = "quasiLikelihoodAdjusted"`. CLI
`singlecell-help` documents its experimental status. The public low-level
`VivoOmicsNBQLInference.fit` also accepts identified gene × pseudobulk counts,
one common design, log offsets, a nonzero contrast and supplied trend dispersions.
It computes every upstream stage itself; cached-fit testing stays internal.

For the one-constraint null `contrast' beta = 0`, the existing constrained NB
solver uses the final scaled trend dispersion. The F statistic is the stable NB
likelihood ratio divided by posterior quasi-dispersion, with numerator DF one.
Denominator DF is adjusted residual DF plus gene-specific prior DF, capped by
the sum of ordinary residual DF over the complete fitted family. Logarithmic
F tails retain probability underflow; an exactly zero LR has P = 1.

The [edgeR manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf)
and pinned edgeR 4.10.5 implementation apply the Poisson bound only to legacy QL.
Modern adjusted-residual QL explicitly disables it. Accordingly, native output
reports `not-applicable-to-modern-adjusted-QL`; this is not a missing modern step.

The result retains original feature identities/statuses, raw unshrunk ML log2
effects, F/P/BH, adjusted residual quasi-dispersion, posterior variance and test
DF. Conditional Wald SE/z/intervals are not presented as QL effect uncertainty;
those expression fields are absent. Conditional fit diagnostics remain available.
`quasiLikelihood.featureIndices` maps local prior/test rows to the original table.
Failed native prerequisites, moderation and null fits retain indexed diagnostics
without manufacturing P values. BH uses only available tests for one contrast.
An optional Cook threshold excludes hypotheses after fitting the complete prior;
it does not change prior membership or the ordinary-DF cap.

The current chain is bounded to 3–100,000 genes and one million gene × pseudobulk
entries, with additional existing design/solver bounds. It never materializes a
dense cell × gene matrix. Active-donor profiles and either effect-prior option
are explicitly rejected with adjusted QL. Legacy QL, varying-support borrowing,
effect intervals, prior uncertainty and selection-adjusted calibration remain open.

## Evidence on experimental data

The frozen [protocol](PROTOCOL.md) uses all 58 original full-support Kang, Hagai
and Crowell arms: 483,576 gene/arm tests. Every native call recomputes the complete
chain and exactly matches the published native refits, adjusted residual DF,
abundance and moderation. No family or gene is dropped. Tight pinned edgeR null
fits and independently recomputed scores, LR, F tails and BH all pass.

| Check | Maximum error | Frozen limit |
| --- | ---: | ---: |
| Constrained means, relative to max(1, reference magnitude) | 3.613e-7 | 2e-6 |
| LR versus tightly converged null, relative | 5.025e-8 | 2e-6 |
| Null scaled score | 9.999684e-8 | 1.001e-7 |
| F arithmetic, relative | 1.122e-15 | 2e-12 |
| Denominator DF / ordinary cap | 0 | 2e-12 |
| Same-statistic log F tail, absolute | 1.225e-12 | 2e-7 |
| BH arithmetic | 6.711e-13 | 2e-12 |

Native CPU FP64 ran on the physical Mac mini. Standalone family runs recorded
0.73–60.14 seconds and maximum resident memory of 140,050,432 bytes. These are
observed run measurements, not controlled comparative speed or million-cell claims.

The existing internal pseudobulk expression owner additionally ran on the
original public Kang aggregation, preserving its counts, design, trend, evidence
and all 15,706 feature identities. All 5,400 eligible tests passed independent
R checks; 906 had BH ≤ 0.05. Runtime was 14.56 seconds and maximum resident
memory was 302,972,928 bytes. The harness links the
actual broad H5AD Omics library with `@testable import` to reach the existing
internal aggregation entry point; it does not introduce an unchecked public API
or constitute a fresh raw-H5AD import. Separate controlled tests exercise the
public dataset route, serialization, failure retention and unsupported options.

All 50 count-model tests in nine suites pass on final source. Both the standalone
inference and broad H5AD library build, as do the full CLI and its help command.
Curated dependency lists were updated for existing scoped build scripts; older
scientific matrices were not rerun merely because their dependency closure grew.

## Calibration remains unresolved

At the original BH ≤ 0.05 threshold, summing across ten sham splits per study:

| Study / dispersion trend | Native adjusted QL | Original default adjusted reference | Original Wald |
| --- | ---: | ---: | ---: |
| Kang / native | 3 | 2 | 0 |
| Kang / edgeR | 3 | 3 | 0 |
| Hagai / native | 39 | 39 | 20 |
| Hagai / edgeR | 79 | 79 | 20 |

These are call inventories on already inspected families, not estimates of FDR
or independent power. Native and original-default references differ upstream;
agreement with independent arithmetic does not establish default equivalence.
Kang sham 05 feature 3809 has native BH 0.03428 versus reference 0.09585. Hagai
treatment feature 9480 has native BH 0.04998 versus reference 0.05004. These are
the only called-set differences across the 58 arms; complete probabilities and
identities are retained. Fresh predeclared calibration and held-out biological
comparisons are required before method selection or default promotion.

## Reproduction and retained attempts

```sh
python prepare.py --abundance-root ABUNDANCE_RUN --moderation-root MODERATION_RUN --ql-root QL_RUN --out RUN/inputs
bash build.sh NATIVE_RUN/build-final
python run.py --inputs RUN/inputs --out RUN/family-final --remote-out NATIVE_RUN/family-final --binary NATIVE_RUN/build-final/nb-ql-inference
python summarize.py --root RUN --ql-root QL_RUN
python archive.py --root RUN --out ARCHIVE
python archive.py --out ARCHIVE --verify
```

`run.py` requires IPv4 SSH to `macmini`, a native executable on that host, and
local R with edgeR 4.10.5 / limma 3.68.5 under the configured library path.
It spools native output remotely, transfers gzip files and checkpoints exact
input/binary/protocol/checker/result hashes. Use a new directory for a failed
attempt; do not overwrite partial results. `build_product.sh` links the H5AD
library; `check_product.py` reconstructs counts from the original source CSR and
calls `product_reference.R` for independent constrained fits and statistics.

The [archive](evidence/2026-09-10/manifest.json) retains all inputs, final outputs,
the initial pilot, build/test/reference logs, source hashes and execution metadata.
It also retains a failed duplicate-source build, the first test failure before
upfront zero-contrast validation, and the first product harness compile failure
before using the required testable import. Final source passed all 50 tests and
all 58 families after the validation repair. Preexisting H5AD `try` warnings are
retained. Archive verification proves artifact integrity; the protocol and
independent checks define the narrower numerical and software qualification.
