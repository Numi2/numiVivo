# File-backed embedding evidence, 2026-09-09

`checks.json` binds the native, frozen numerical-oracle and independent reference
results through exact result and graph hashes. `runtime.json` records the source
base, executables, hosts and measurement boundaries. `production-sources.json`
identifies the tested source tree; `protocol.json` predates the real runs.

- `real/`: full-cohort result, graph/PCA state, receipts, native timings and frozen
  numerical-oracle output. Original public H5ADs remain external to this archive;
  their retained source metadata and earlier benchmark evidence identify them.
- `graph-reference/`: every binary coordinate and FP64 value, PCA score byte and
  cell identity compared with the previous independently qualified full graph.
- `reference/`: full schedule counts, curve checks, all three umap-learn optimizer
  coordinate sets, exact quality-query rows and descriptive quality measurements.
- `fixtures/`, `fixture-inputs/`, `graph-fixtures/`, `clustering-regression/`:
  lifecycle, fitted/query, exact/HNSW, 2D/3D, replay and expected rejection evidence.

Large JSON, binary and log artifacts are deterministic gzip (mtime zero).
Decompress them before calling native verifiers or reference scripts. Rebuilding
requires original H5ADs and fresh PCA/graph receipts for the chosen binary;
never relabel an older executable's receipts. Raw compiler logs preserve their
original whitespace.

The frozen oracle preserves the old numerical algorithm and result shape but
uses the current public options validation, including larger explicit budgets.
It does not prove that the earlier CLI accepted those enlarged-budget plans.
Norman uses all cells for optimization and only the predeclared 2,048 queries for
rank-based quality; each query is compared with the whole cohort. This is not
million-cell, convergence, biological or GPU/scverse performance qualification.

`retired-inputs-complete.json` records redundant input snapshots relocated to
verified local APFS clone archives to free remote disk space. Canonical public
sources and previous complete clustering bundles stayed on the Mac mini. Restore
those retired paths from their exact local archives before replaying older bundles.
No original dataset was deleted.

See [file embedding](../../FILE_EMBEDDING.md) for commands, measurements and limits.
