# File-backed clustering evidence, 2026-09-09

`checks.json` aggregates the actual native, frozen-oracle and independent
reference results. `runtime.json` records the base revision, executing binaries,
hosts, measurement scope and retained qualification failures. The production
source manifest identifies every source file compiled into the tested tree.
`protocol.json` was written before the full-cohort clustering runs.

Norman's independent reference job ran during native replay and omitted the
optional frozen-oracle argument. Its corresponding field is explicitly null.
The completed native runner separately proves byte equality; the aggregate
archive check binds that result and the independent checks to the same result
and graph hashes. Raw compiler logs preserve their original whitespace.

- `real/`: per-cohort clustering result, parent graph/PCA state, receipts, native
  lifecycle logs and frozen pre-refactor result. The large original public H5AD
  files are omitted; retained PCA source metadata and earlier dataset evidence
  identify them. Reproduction must refit those original sources with the chosen
  binary. An old receipt cannot be relabeled for a newly built executable.
- `graph-reference/`: every binary coordinate and FP64 value, PCA score byte and
  cell identity compared to the earlier independently qualified full graph.
- `reference/`: independent objectives/connectivity, all three NetworkX runs,
  partitions, reference-to-reference stability and descriptive metadata metrics.
- `fixtures/`, `graph-fixtures/`, `fixture-inputs/`: fitted/query, exact/HNSW,
  sparse format, deterministic replay and expected rejection evidence.

Large JSON, binary and log artifacts use deterministic gzip (mtime zero).
Decompress them before passing bundles to native verification or checkers.
The original full-cohort H5AD inputs are intentionally external to this archive.
The three full native workflows were sequential on the Mac mini; references
ran on the laptop. Native publish and verify measurements include repeated PCA
and graph reconstruction. Frozen-oracle timing has a different workload and
must not be treated as a workflow speed comparison.

See [file clustering](../../FILE_CLUSTERING.md) for measured values, commands,
resource bounds and remaining work. This is numerical storage/algorithm
qualification, not biological, million-cell, Leiden, GPU or scverse performance
qualification.
