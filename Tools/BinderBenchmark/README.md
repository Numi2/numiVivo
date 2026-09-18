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

Portable regression (Swift 6, Python 3 standard library):

```sh
python3 Tools/BinderBenchmark/check_native.py
```

This compiles the exact owner and its XCTest source in a temporary package; it does
not build or qualify the full Apple/Metal product. On macOS the same tests are part
of `swift test --filter VivoBinderBenchmarkTests`.

Next: pinned experimental-table import, native CLI, structure-derived observations,
then independently qualified preparation/MD. No NVIDIA or paid service is required
for this first retrospective scoring step.
