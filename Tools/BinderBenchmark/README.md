# Binder selection: native retrospective benchmark

`VivoBinderBenchmark` in `NumiVivoKit/Binder` owns binary-outcome candidate evaluation.
It is not a molecular simulator, affinity predictor or validated binder-selection product.

A dataset names one assay endpoint, exact source SHA-256, grouping method and each
candidate's target, group, raw outcome, preserved source fields and available numeric
features. Unknown/untested/inconclusive/expression-failure outcomes remain nonbinary.
The caller must preserve and verify the source bytes; a supplied digest is not proof
that the data was faithfully imported.

An explicit plan separates whole training targets from whole test targets. Training
candidates whose declared sequence groups occur anywhere in the test set are purged.
Only training rows determine feature scaling and the ridge-logistic model. Model and
baseline use the same complete-case candidates; every exclusion remains in the report.
The fixed baseline feature is a higher-is-better score, not a probability.

Reports include fitted parameters/convergence, exact candidate IDs, held-out predictions,
AUROC, threshold-block average precision, Brier score where meaningful and expected
hits/precision at K. Ties receive fractional selection credit, independent of row order.
Single-class AUROC remains unavailable. Reports do not pool targets into a potentially
misleading aggregate or grant a qualification based on a positive-looking result.

Freeze the plan before inspecting test outcomes. A code path cannot establish that
historical fact. Exact-sequence grouping does not establish homology independence;
related targets/designs are not independent experimental replicates. Start with a
score-only baseline/ensemble before testing any additional structural/MD feature.

Portable regression (Swift 6, Python 3.10+ standard library; OpenSSL development library on Linux):

```sh
python3 Tools/BinderBenchmark/check_native.py
python3 Tools/BinderBenchmark/test_tools.py -v
```

This compiles the exact owner and its XCTest source in a temporary package; it does
not build or qualify the full Apple/Metal product. On macOS the same tests are part
of `swift test --filter VivoBinder`.

Current: published-table import, native CLI and source-replayable score evaluation.
Next: structure-derived observations, then independently qualified preparation/MD.
No NVIDIA or paid service is required
for this first retrospective scoring step.

## Published-source workflow

The adapter targets Anthropic's `claude-protein-binder-design` release (CC BY 4.0),
revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`,
`data/tables/design_summary.csv`, Git blob `d1573ba03e8322c70ccb3a40e46418e86b40e2dc`.
Source: https://huggingface.co/datasets/Anthropic/claude-protein-binder-design

The downloader checks that upstream blob identity before publishing. The importer
itself verifies input consistency and fingerprints; it does not authenticate an
arbitrary local CSV's upstream origin. Preserve the downloader's `SOURCE.json`.
Initial scope is BBF-14, MBP and EGFR, explicitly selected rather than all targets.
Each vendor is a separate benchmark. Only exact `binder` and `non_binder` source
classes become binary; other classifications remain unavailable with raw text.
No combined `binder_final`, affinity, expression or assay-derived field enters the
predictors. All source fields remain available for review.

```sh
python3 Tools/BinderBenchmark/fetch_source.py /tmp/numi-binder-source
swift build -c release
BIN=.build/release/numivivo
$BIN binder-import /tmp/numi-binder-source/source.csv \
  /tmp/numi-binder-source/import-adaptyv.json /tmp/numi-binder-input
python3 Tools/BinderBenchmark/make_plan.py /tmp/numi-binder-input \
  --test-target BBF-14 --output /tmp/numi-binder-plan.json
$BIN binder-evaluate /tmp/numi-binder-input /tmp/numi-binder-plan.json /tmp/numi-binder-result
$BIN binder-verify /tmp/numi-binder-result
```

The example compares a preselected Boltz-2 ipSAE ranking with a training-only
three-score ensemble and training-prevalence baseline. This is not Numi physics
improving binding prediction. It establishes the benchmark against which added
physical observations must later be tested. Do not tune the plan on the test result.

Bundles retain original CSV/config bytes, deterministic imported records, a plan
and complete report where applicable. SHA-256 binds every file and the producing
executable. Verification reconstructs the import and repeats evaluation; modified
records cannot pass merely by updating a file checksum. Rebuilt executables have a
new identity and require a fresh run, not relabeling old evidence. Outputs must be
new directories; incomplete runs do not publish and existing data is never replaced.

The portable checker uses the existing `Tools/Posterior/PortableSupport.swift`
OpenSSL/canonical-JSON facade. It tests the exact binder source and CLI, not the
full production artifact store, Apple package integration or Metal execution.

## Executed evidence and next boundary

The [2026-09-18 portable record](evidence/2026-09-18-portable/README.md) records
35 passing Swift tests, seven passing offline tool tests and a 90-candidate synthetic
CLI campaign with independent Python checks. This is implementation evidence only.
The full published CSV was not downloaded/scored in that environment, and the
complete macOS/Metal application was not built. The actual source-file downloader
was exercised with mocked bytes, not the remote endpoint.

This increment establishes a score-only comparison, not the previously proposed
AI-structure-to-MD pipeline. It does not yet provide structure repair, protonation,
parameter assignment, solvent setup, physical interface descriptors, inference from
Boltz-2 weights or improved selection on measured outcomes. These should be added
through the existing molecular owners and measured against the fixed baseline,
not reimplemented as another simulation subsystem.
