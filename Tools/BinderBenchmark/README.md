# Binder selection: native retrospective benchmark

`VivoBinderBenchmark` in `NumiVivoKit/Binder` owns binary-outcome candidate evaluation.
It is not an affinity predictor or a validated binder-selection product.

**Current execution:** the pinned published table has been downloaded and the fixed
six-fold score-only campaign completed, with independent numerical checks and a
separate local rebuild reproducing every complete evaluation report. The three-score
ensemble ties the fixed Boltz-2 ranking in four top-10 comparisons and loses in two;
it improves none. This is a negative development result, not evidence of Numi
physical-analysis improvement. See the [published record](evidence/2026-09-18-published/README.md).

**Current implementation:** published-table import, outcome-blind fixed rankings,
explicit training-support checks, raw PDB/mmCIF reconstruction and
[structure-derived geometry with matched score-only controls](STRUCTURAL_ANALYSIS.md).
The new [per-prediction surface analysis](PREDICTION_ANALYSIS.md) computes source-bound
burial and radius-overlap observations. Its fixed 18-structure numerical panel passed
independent reference checks; this is not the full 270-design structural-selection
experiment. The preparation/MD route and full Apple/Metal application remain separate.
See [ranking/support evidence](evidence/2026-09-18-ranking-support/README.md) and
[prediction-surface evidence](evidence/2026-09-18-prediction-surface/README.md).

## Evaluation contract

A dataset names one assay endpoint, exact source SHA-256, grouping method and each
candidate's target, group, raw outcome, preserved source fields and available numeric
features. Unknown, untested, inconclusive and expression-failure outcomes remain
nonbinary. A supplied digest alone does not authenticate a local dataset's origin.

An explicit plan separates whole training targets from whole test targets. Training
candidates whose declared sequence groups occur anywhere in the test set are purged.
Only training rows determine scaling and the ridge-logistic model. Methods use the
same complete-case candidates; exclusions remain in reports. The fixed baseline
feature is a higher-is-better score, not a probability.

Reports retain fitted parameters/convergence, exact candidate IDs, held-out predictions,
AUROC, threshold-block average precision, Brier score where meaningful and expected
hits/precision at K. Ties receive fractional selection credit independent of row order.
Single-class AUROC stays unavailable. No favorable metric automatically grants
qualification, and target/assay endpoints are not pooled into an apparent replication.

Freeze plans before evaluation. Code can bind a plan to results but cannot prove that
nobody previously inspected outcomes. Exact-sequence grouping does not establish
homology independence. The already examined published targets are development data
for subsequent feature work, not a fresh external test.

## Published-source workflow

The adapter targets Anthropic's `claude-protein-binder-design` release (CC BY 4.0),
revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`,
`data/tables/design_summary.csv`, Git blob `d1573ba03e8322c70ccb3a40e46418e86b40e2dc`.
Source: https://huggingface.co/datasets/Anthropic/claude-protein-binder-design

The downloader checks the pinned Git-blob identity before publishing. The importer
checks consistency and fingerprints, not upstream authenticity of arbitrary local
CSV bytes. Preserve the downloader's `SOURCE.json` and original source.

The fixed campaign explicitly selects BBF-14, MBP and EGFR: 270 designs from the
1,440-row source. Each vendor remains a separate benchmark. Only exact `binder`
and `non_binder` classifications become binary. Other text remains preserved and
excluded; neither combined `binder_final`, measured affinity, expression nor other
assay-derived columns are admitted as predictors.

Full fixed campaign using the normal Apple executable:

```sh
python3 Tools/BinderBenchmark/fetch_source.py /tmp/numi-binder-source
swift build -c release
python3 Tools/BinderBenchmark/run_published.py /tmp/numi-binder-source \
  .build/release/numivivo /tmp/numi-binder-campaign
```

All output directories must be new. The campaign writes all six plans before
import/evaluation, then checks source identity, reconstructed import, training
population/scaling/fit gradients, held-out predictions and metrics independently
in Python. It preserves the original CSV, executable, plans, bundles, logs and hashes.

Individual native commands:

```text
numivivo binder-import SOURCE.csv IMPORT.json NEW_BUNDLE
numivivo binder-ranking-query IMPORT_BUNDLE NEW_QUERY
numivivo binder-rank QUERY_BUNDLE RANKING_PLAN.json NEW_RANKING
numivivo binder-assess-ranking IMPORT_BUNDLE RANKING_BUNDLE NEW_ASSESSMENT
numivivo binder-evaluate IMPORT_BUNDLE PLAN.json NEW_RESULT
numivivo binder-evaluate-structures IMPORT_BUNDLE PLAN.json STRUCTURES.json NEW_RESULT
numivivo binder-evaluate-supported IMPORT_BUNDLE PLAN.json POLICY.json NEW_RESULT
numivivo binder-evaluate-structure-sources IMPORT_BUNDLE PLAN.json SOURCES.json NEW_RESULT
numivivo binder-analyze-prediction PREDICTION.json NEW_BUNDLE
numivivo binder-verify BUNDLE
```

Use `make_plan.py IMPORT_BUNDLE --test-target TARGET --output NEW_PLAN.json` for a
score-only plan. Structural evaluation requires explicitly adding selected geometry
features to a separately frozen plan; see [the structural contract](STRUCTURAL_ANALYSIS.md).

## Reproducibility and tests

Bundles retain source/configuration bytes, deterministic imports and complete reports.
SHA-256 binds files and the producing executable. Verification reconstructs import,
geometry where applicable, fitting and evaluation. Replacing an artifact and merely
updating its recorded checksum cannot make inconsistent results pass. Rebuilt
executables have new identities and need fresh runs, not relabeled old evidence.
Existing outputs are not overwritten and incomplete runs do not publish bundles.

```sh
python3 Tools/BinderBenchmark/check_native.py
python3 Tools/BinderBenchmark/test_tools.py -v
```

Requirements: Swift 6, Python 3.10+ standard library and, on Linux, the OpenSSL
development library. `native_package.py` compiles the exact binder, molecular
structure, validator, interface and artifact-primitive sources with the native CLI.
Only canonical JSON and hashing use `PortableJSON.swift`, a harness-only
Foundation/CryptoKit-or-OpenSSL implementation. This is not the full production
artifact store or Apple package. On macOS the native tests are also included in
`swift test --filter VivoBinder`.

A retained portable executable can run the same published campaign:

```sh
python3 Tools/BinderBenchmark/native_package.py /tmp/numi-binder-native --test
python3 Tools/BinderBenchmark/run_published.py /tmp/numi-binder-source \
  /tmp/numi-binder-native/.build/debug/BinderCLI /tmp/numi-binder-portable-campaign
```

The focused Binder benchmark CI at `cbcf613b5a283f2a45b28ef010006195f5b0e984`
passed 50 Swift tests, seven Python tests, both synthetic CLI campaigns and the real
six-fold campaign. This does not assert success of unrelated repository workflows.
The [earlier portable record](evidence/2026-09-18-portable/README.md) remains a
historical account of the preceding 35-test, synthetic-only implementation.

## Remaining scientific and engineering work

Source-bound per-prediction analysis is implemented and numerically checked on
18 real structures. Next are complete cohort processing, real multi-model/seed
confidence and pose observations, and a frozen structural-feature comparison against
both fixed rankings. The compact SDK accumulator does not yet provide a persistent
series CLI or independently evaluated learned structural correction. Retained source
and sequence identity do not authenticate unknown model execution details.

Structure repair, protonation, parameter assignment, solvent setup, qualified
molecular sampling, model-weight inference and measured Apple acceleration are not
provided by this increment. Integrate those through existing molecular owners;
do not introduce another simulation subsystem or treat static geometry as an energy.
