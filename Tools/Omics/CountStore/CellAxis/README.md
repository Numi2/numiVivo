# File-backed cell identities, QC and count aggregation

The native cell axis stores cell identities and per-row declarations on disk.
The count consumer streams sparse records, writes cell QC and aggregate membership
to disk, and retains only bounded feature/sample/group dictionaries and aggregate
counts in memory. These are product CLI commands, with Swift owners in
[`VivoFileCellAxis.swift`](../../../../Sources/NumiVivoKit/Omics/VivoFileCellAxis.swift)
and [`VivoFileCountStream.swift`](../../../../Sources/NumiVivoKit/Omics/VivoFileCountStream.swift).

## Evidence as of 2026-09-11

| Check | Result | Scope |
|---|---|---|
| Native regression tests | PASS: 14 tests in 2 suites | Ten new file-axis tests and four existing stream tests |
| Complete Parse cell-axis import | PASS: 725,031 cells | Every original barcode/sample byte, row cardinality and exact retained-matrix total |
| Native axis reopen | PASS | Byte-identical receipt after full structural validation |
| Native import peak RSS | 65,617,920 bytes (62.6 MiB) | Native child process; excludes Python adapter and OS file cache |
| Import elapsed time | 13.01 seconds | Includes input streaming waits; not isolated kernel throughput |
| Full paired count ingestion and source replay | PASS: both complete phases | Same 1,373,870,697 records; exact QC, membership and aggregate agreement; [full results](COUNT_RESULTS.md) |
| Biological prediction | Not fitted or scored | Storage qualification does not establish an outcome prediction |

The [retained cell-axis evidence](evidence/2026-09-11-axis/summary.json) includes
the native artifact, executable, 105 build-input identities, tests, source
bindings and exact executed preparation recipe. All 6,916 source bindings were
checked against the previously retained Parse preparation and completed count
archives. Declared totals come from the retained matrix, whose counts differ from
historical source QC; no historical totals were substituted.
Fresh extraction of all 23 archived files and
[native reopening of the restored axis also pass](evidence/2026-09-11-restoration/summary.json),
with every restored hash and the original receipt reproduced exactly. This checks
artifact recovery, not another count-stream replay.

The full count comparison is a separate gate. It requires identical complete
canonical byte streams, all aggregate coordinates, every cell QC record and
every aggregate membership to agree with the previous native implementation,
fresh independent sums and the earlier frozen source qualification. Both
ingestion and a new source replay must pass before reporting its memory comparison.
The earlier donor-partitioned process measurements cannot substitute for this
same-cohort comparison. Frozen pending count recipes are retained explicitly;
active or partial count outputs are excluded from this axis archive.

## Use

```sh
numivivo singlecell-cell-axis-import --header header.json --output axis < cells.jsonl
numivivo singlecell-cell-axis-verify axis
numivivo singlecell-file-count-stream --axis axis --output counts < records.bin
numivivo singlecell-file-count-stream-verify counts < records.bin
```

The header contains schema version 1, `metadata` with the existing feature and
sample dictionaries and an **empty** `cells` array, `cellCount`, `nonzeros`,
`hasRowTotals`, an `annotations` dictionary and `sourceDeclaration`.
Each newline-terminated JSONL row has `barcode`, `sampleID`, `nonzeros`, optional
`group`, and `totalCounts` when declared by the header. Unknown top-level header
and row keys are rejected. Source decoding and provenance remain the adapter's
responsibility; these commands do not themselves parse arbitrary HDF5 metadata.

Count records use the existing 16-byte little-endian format: unsigned 32-bit
global cell index, unsigned 32-bit feature index, unsigned 64-bit positive count.
Rows/features must be strictly ordered without duplicate coordinates. Exact
declared cardinalities and optional totals are checked, including empty cells;
truncation, invalid indices and integer overflow fail without publishing a bundle.

## Disk format and memory bounds

`axis/rows.bin` holds 40 bytes per cell: string offset (u64), barcode byte length
(u32), original sample-reference byte length (u32), original group byte length
(u32; zero means absent), sample index (u32), annotation index (u32; `UInt32.max`
means absent), nonzero count (u32), and total counts (u64).
`axis/strings.bin` concatenates the original UTF-8 barcode, sample reference and
optional group text. Canonically equivalent Unicode spellings remain byte-exact;
duplicate sample/barcode identities are still rejected using Swift string identity.
`UInt64.max` is a valid total, not a missing-value sentinel.

The count bundle embeds the axis and adds `quality.bin`, `report.json` and
`receipt.json`. Each 24-byte QC row contains total counts (u64), mitochondrial
counts (u64), nonzeros (u32) and group ordinal (u32). Thus all cell memberships
remain addressable without resident membership arrays. Groups report the explicit
`sourceCellCount`; aggregation preserves technical-sample pooling, annotations,
empty groups and the existing biological-replicate/condition semantics.

Headers are limited to 64 MiB, identifiers to 1,024 UTF-8 bytes, JSONL rows to
16 KiB, cells to 20 million, source records to 40 billion, and features to
100,000. Samples and annotation dictionaries are each limited to 100,000 entries.
The count consumer caps observed groups at 100,000, sample/group associations
at one million and aggregate nonzeros at five million. These are admission bounds,
not demonstrated performance at every maximum. Reads and output buffers are
bounded; transient duplicate detection uses a disk hash table. Private snapshots
use APFS clones where available and bounded copying otherwise.

Receipts bind header, row, string, QC, report and canonical input bytes. Opening
a snapshot checks immutable bytes and structure; full numerical reconstruction
requires the explicit count-stream verification command with the original stream.
A valid hash alone does not prove numerical correctness.

This removes the resident cell/QC/membership requirement from this count path.
Other H5AD, PCA, graph and prediction paths still have their existing bounds.
Whole-pipeline out-of-core execution and Metal acceleration are not established
by this qualification.

## Reproduce and retain

Build the scoped CLI with [`../../H5AD/build.sh`](../../H5AD/build.sh), then run
[`test.sh`](test.sh) with its build directory. The Parse tools are deliberately
bound to the dated, already-qualified source preparation rather than a generic
dataset convention:

- [`prepare_parse.py`](prepare_parse.py) streams the complete original cell axis
  and exact matrix totals into the native import, then independently checks all rows.
- [`run_parse_counts.py`](run_parse_counts.py) tees each canonical record to both
  native consumers and records each child's process metrics. It checks fetched
  range identities against the earlier completed count qualification.
- [`compare_parse_counts.py`](compare_parse_counts.py) compares complete native
  reports, all QC rows, memberships and independent matrices with bounded readers.
- [`retain_axis.py`](retain_axis.py) packs completed cell-axis evidence only.

The archive uses the existing content-addressed `members.json`/`objects/<sha256>`
format. Verify or restore with `Tools/Omics/PerturbationPrediction/Duration/archive.py`.
See [dependencies.json](evidence/2026-09-11-axis/dependencies.json) for the exact
source archives needed to repeat independent source checks. Restoring bytes is
not another native execution or network replay. Parse data and derivatives are
subject to **CC BY-NC 4.0**.

## Complete count and analysis handoffs

[Both full same-cohort count phases now pass](COUNT_RESULTS.md): 1.37 billion records, exact cell/QC/membership/aggregate agreement, and native peak RSS of 165.6 MiB for ingestion and 181.8 MiB for source replay. Count execution, offline archive checks and restored-artifact checks are separate evidence stages.

The [file-backed expression owner](../Expression/README.md) now consumes these aggregates through explicit receipt-bound membership references and reproduces the prior native statistics on the complete cohort. All six external statistical comparisons complete. This closes the interface without fabricated cell lists; the full analysis report still peaks near 1 GiB and no Parse prediction is fitted or scored.
