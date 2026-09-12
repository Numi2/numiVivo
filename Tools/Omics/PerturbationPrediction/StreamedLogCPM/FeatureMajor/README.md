# Feature-major streaming normalization

The log-CPM plan accepts `"ordering":"featureMajor"` for strictly ordered
(feature, row) positive count records. Omitted ordering remains `rowMajor`.
Unknown orderings are rejected. Each cell's observed UInt64 count total is
accumulated with overflow checks and must equal its declared full-axis total
before a result is returned. Duplicate or decreasing coordinates, out-of-axis
records, excess/missing counts and reuse after failure are rejected.

This consumes the existing Norman source-major store directly, without sorting,
a second multi-gigabyte count copy or a dense cell-by-gene matrix. It retains
row totals, assignments, and condition-by-feature sums/corrections in memory.

All 111,445 cells, 33,694 features, 237 conditions and 361,582,621 records passed.
Every one of 7,985,478 normalized means matches independent NumPy indexed
accumulation within 9.06e-14 absolute error. The reference also rechecked the
complete count digest, strict feature-major ordering and all row totals.
The actual scoped native CLI exactly reproduces the standalone result and
rejects incorrect/unknown ordering declarations. Small controls verify layout
equality, zero-cell denominators, malformed input and closed-state behavior.

The standalone full-data run took 5.55 seconds and peaked at 766,869,504 bytes
RSS. Large condition-by-feature result arrays and JSON serialization remain
resident; this is not million-cell or fully disk-backed-output qualification.
No biological prediction or independent replicate claim follows from these
normalization checks. Norman pooled perturbation conditions are retained as
conditions, not relabeled as independent biological donors.

[verification.json](verification.json), [cli-receipt.json](cli-receipt.json),
[protocol.json](protocol.json), [run.log](run.log) and [cli.log](cli.log) retain
source-bound results and process measurements. [reference.py](reference.py),
[check_cli.py](check_cli.py) and [Errors.swift](Errors.swift) retain the executed
verification routines with original host paths.

Complete native outputs, the reference array, plan and binaries remain on the
execution host at `/Users/n/numivivo-logcpm-feature-major-20260912`. Their hashes
are bound in these receipts. They are not embedded in this documentation bundle.
Reproduction requires the separately retained full Norman count store and
metadata at `/Users/n/numivivo-count-store-final-20260909/store`, or reconstructing
that [qualified store](../../../CountStore/README.md). Use a new output directory;
the checks refuse to overwrite previous completed outputs.
