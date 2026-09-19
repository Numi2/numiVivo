# Prediction surface numerical qualification — 2026-09-18

Execution occurred September 18 in America/Los_Angeles (September 19 UTC).
This is structural numerical qualification, not improved experimental selection.

## Delivered

`24114a44` extends the existing molecular-interface owner with explicit-radius
heavy-atom surface burial and packing observations. `e5ebe631` adds outcome-blind
per-prediction analysis, source-replayable CLI bundles and compact SDK sample
summaries. `015767cd` strengthens cross-predictor identity checks, adds the real
published-data checker and its dedicated CI workflow, and tests quadrature error.

No separate molecular parser, MD engine or artifact store was introduced.
The [usage and scope](../../PREDICTION_ANALYSIS.md) describe the supported fields,
units, sampling limits and remaining batch/PAE/ranking integrations.

## Real source processing

The existing hash-selected panel supplied six designs from each of BBF-14, MBP and
EGFR, with one seed-best Boltz-2 1:1 prediction per design. All 18 source predictions
were available and passed, covering 62,660 heavy atoms. The original 5,534,313 CIF
bytes, construct FASTA, manifest and source table were retained and verified.
Selection and chain matching did not use experimental outcomes.

The published source is Anthropic `claude-protein-binder-design`, CC BY 4.0,
revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`. Acquisition artifact:
https://github.com/Numi2/numiVivo/actions/runs/35415237915/artifacts/10575506988
ZIP SHA-256: `b5911cbb989cf818c9f4e6c7b8bbc724c9ca695b4fcda3dd1aba028f68f504ee`.

At 960 sphere points per atom with the explicit CNOS/probe profile in
[protocol.json](protocol.json), native areas were independently checked against
Biopython 1.86, and radius-overlap counts against NumPy 2.3.5. All 18 comparisons
and subsequent native source-replay verifications passed.

Maximum observed per-atom area difference was 0.0012579460583749258 nm^2 (within
one quadrature point). Maximum total summed-burial difference was
5.861977570020827e-13 nm^2. All overlap pair counts agreed. The comparison allows
explicit atomic/total tolerances because reference coordinates/sphere points and
native arithmetic have different precision. These observations do not establish
convergence to a continuous surface or experimental binding correctness.

The compact [summary](summary.json) retains every candidate, source hash, result and
numerical difference. Reported local publication times include the native
calculation and its automatic publication replay; they are not inference timings
or a cross-hardware performance benchmark.

## Executed regression scope

The exact native source subset through `015767cd` passed 109 Swift tests, zero
failures. Thirty-six offline Python tests passed: seven original source/plan tests,
ten source-transport tests, ten panel-acquisition tests and nine new checker tests.
The full native harness also passed the existing 90-row score campaign, the
17-structure embedded and raw PDB/mmCIF campaign, independent geometry/fit checks
and all mocked-source folds. The outcome-blind ranking CLI checker passed.

All six real legacy fitted evaluations were rerun with the new release executable.
Their complete report bytes remained identical to the previous ranking-support
replay archive, not merely their rounded metrics. No baseline was retuned.

Local Swift 6.2.1 ran on x86_64 Linux, Python 3.13.5. The native owners and CLI were
compiled with the existing harness-only canonical-JSON/OpenSSL facade. This is not
a build of the complete Apple package, production artifact store or Metal backend.
Source/test hashes and executable identity are retained in [execution.json](execution.json).

The existing binder CI completed successfully at `015767cd` in run
https://github.com/Numi2/numiVivo/actions/runs/35417677649 (job `105829349318`).
Its published baseline, fixed-ranking and support-decision steps all passed.
The separate surface-reference workflow also completed successfully in run
https://github.com/Numi2/numiVivo/actions/runs/35417677604 (job `105829349320`).
Its downloaded artifact `10575849580` passed archive-digest verification. All 18
complete native analysis reports matched the local reports byte-for-byte, despite
distinct producing executable identities. Both runs independently compared original
source structures with the reference implementation. This is agreement for these
runs, not a general cross-platform bitwise-equivalence claim.

Surface CI artifact ZIP SHA-256:
`1bdbce0850495d6e8ed8bcb2bb623b11ec37463b3aa71880b9fc587ae52aeaff`.
The separately downloaded source snapshot from the successful binder CI matched
all 18 compiled source files, 11 native test files and the Python checker files.
Both CI scopes remain subsets of the product, not a full Apple or Metal build.

## Scientific boundary

Actual evidence is 18 seed-best Boltz-2 predictions. Multiple predictors/seeds and
contact-agreement summaries were exercised in synthetic SDK tests, not a complete
real multimodel campaign. The 270-design structural comparison, PAE/token mapping,
pose alignment, persistent series CLI, learned structural correction and independent
experimental evaluation remain unfinished. Surface features are not yet inputs to
a newly validated candidate selector.

No improved top-ten selection, inferred binding affinity, hydrogen-bond/polarity
assignment, model-weight inference, MD result or Apple speedup is claimed.
