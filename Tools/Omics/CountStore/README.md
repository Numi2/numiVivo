# Persistent count storage and windowed normalization

`singlecell-h5ad-store` retains an exact source snapshot, explicit mapping,
count metadata and row/feature statistics alongside a binary sparse count file.
Counts are written in 1 MiB payload buffers. Normalization consumes immutable
private snapshots through 16 MiB POSIX mapping windows and unmaps each previous
window. Count payloads are not materialized as Swift arrays or JSON matrices.
Metadata, library totals and feature statistics remain resident, as does one
sparse source segment during AnnData canonicalization.

## Product commands

```sh
numivivo singlecell-h5ad-store original.h5ad --plan mapping.json --output raw-store
numivivo singlecell-count-store-verify raw-store
numivivo singlecell-count-store-normalize raw-store --target 10000 --output normalized
numivivo singlecell-count-store-normalize-verify normalized --store raw-store
```

Import uses the existing AnnData count decoder, including explicit count-layer
selection, CSR/CSC/dense admission, exact integer validation, duplicate-coordinate
summation, zero removal and identity mapping. The stored traversal order follows
the source's major axis. Two source layouts can represent the same matrix with
different record order; a repeat of the same source is deterministic.

Each raw record is exactly 16 bytes: little-endian UInt32 row, UInt32 feature,
UInt64 count. Normalized records have the same coordinate fields followed by
Float64 `log1p(count * target / libraryTotal)` instead of UInt64 counts. Their
format identifiers distinguish those meanings. The raw values are never replaced
by floating-point approximations. Source and normalized receipts bind metadata,
record payloads, parameters and implementation separately.

This first store path covers the existing H5AD count metadata interface. It is
not yet a disk-backed replacement for every multimodal or downstream algorithm.
The record window primitive can support later native consumers; it is not a new
scientific inference engine or a second source of count semantics.

## Reconstruction and bounded I/O

Publication requires a new destination and removes private staging files on
failure. Snapshot copying uses a fixed POSIX buffer and handles partial writes
and interrupted calls. File descriptors reject symlinks and non-regular files.
Mapped files belong to private snapshots so an external truncation cannot invalidate
an active window. File sizes and entry bounds are checked before mapping.

Raw verification reconstructs metadata, statistics and the complete record digest
from retained source and mapping. Normalization first performs that verification,
then reads the count store through mapping windows. Normalized verification
recomputes the transform and checks all receipt-bound artifacts. Command timings
therefore include source reconstruction. Rehashing modified records cannot bypass
it. Verification hashing does not require a second resident count matrix.

Current engineering limits are two million observations, 200,000 features,
two billion entries, 64 GiB per source/snapshot file, 512 MiB encoded metadata and
256 MiB encoded statistics. These are admission limits, not qualifications at all
of those sizes. A 1 MiB record flush temporarily has a Data copy in addition to
its write buffer. Window size does not describe total process or OS cache memory.

## Qualification protocol

`check_interop.py` consumes the existing H5AD interoperability fixtures. It checks
all CSR/CSC/dense entries and integer values above 2^53, exact marginal statistics,
FP64 normalization, source reconstruction, repeatability, invalid inputs and
rehashed count/normalized-value tampering. Swift controls cross a mapping-window
boundary, revisit an earlier window and test empty/truncated/non-regular inputs.

`run_native.py` operates on the complete measured Norman filtered release on the
native build host, retaining all 111,445 cells and 33,694 genes. The source SHA-256
is `efde6f5301fe256725dce1d980f37bd96a13481a9a16135515897368e631affc`.
The source and mapping provenance are documented in the [Norman benchmark](../PerturbationPrediction/Norman/README.md).
The native benchmark gates each command below 1 GiB maximum resident memory on
that dataset. `check_norman.py` compares every remote record with local SciPy and
Scanpy, streaming large files over SSH rather than making another local disk copy.
Different reference/native hosts preclude a comparative speed claim.

The first full run completed reconstruction but normalization peaked at 6.79 GB
and normalized verification at 12.58 GB. This footprint tracked the snapshot
files copied through the Foundation read path. Replacing it with a fixed POSIX
buffer reduced normalization to 282 MB on the same complete dataset. Both the
initial measurements and subsequent qualification remain in evidence. There is no million-cell,
parallel-kernel, downstream graph/integration or Metal qualification from this
store milestone.

## Complete Norman native result

The 2026-09-09 release on Apple M4 Pro / macOS 26.6 retains 111,445 cells,
33,694 genes and 361,582,621 nonzero records. Each count/value payload is
5,785,321,936 bytes. Native source reconstruction passed with the following
end-to-end measurements (decimal MB):

| Command | Wall seconds | Maximum resident MB |
| --- | ---: | ---: |
| Import | 41.24 | 254.75 |
| Normalize, including source verification | 64.35 | 281.92 |
| Verify normalized store, including source and transform reconstruction | 70.09 | 285.07 |

The selected Swift gate passed 65 tests in 18 suites. The final CLI interoperability
gate passed 26 commands including 13 expected rejections. These include structural
fixtures; the full Norman release supplies the measured-data scale check.

[Retained evidence](evidence/2026-09-09/) includes the initial failed memory gate,
final tests, native command logs, receipt hashes and exact source/product bindings.
The qualified product SHA-256 is
`ab43b10b5b4761afe0d56769f3f794e8a074a51c030f9564d2196d27cd6aa70b`.

The independent checker compared every raw and normalized record. All counts,
coordinates, cell/gene identities and marginal statistics match SciPy exactly.
All 361,582,621 normalized values match Scanpy FP64 normalization with maximum
absolute error `8.881784197001252e-16` (fixed absolute/relative tolerance `1e-12`).
[The final receipt](evidence/2026-09-09/final/checks.json) includes both complete
payload digests and the reference environment. Lossless gzip transport reduced
SSH transfer volume; every decompressed record was checked and hashed.

## Explicit Metal normalization and retained historical payloads

The [Metal qualification](Metal/README.md) adds an explicit `metal-fp32` backend,
with exact coordinates, a declared FP32 numerical profile and device-bound
reconstruction. The complete original Kang cohort passes independent checks
and native replay. Its three-run end-to-end median does not beat CPU; the
historical FP64 default and complete CPU output bytes remain unchanged.

Both historical Norman `normalized/values.bin` payloads (initial and final runs)
were identical. They now share one verified local gzip backup, with every decoded
byte checked against both receipts before deleting the redundant remote files.
Their original receipts, metadata, logs and memory-failure evidence remain.
The [restoration manifest](Metal/evidence/2026-09-11/manifest.json) records exact
paths and compressed/decoded hashes. Restore the values file from that backup
before replaying either historical normalized bundle. Later benchmark duplicate
payloads have separate retained-copy mappings in the same archive.
