# Storage access at full HIRISA dimensions

Large shuffled CSR reads now use a 4 KiB positional-read buffer. Files that fit
within one 16 MiB window keep the existing mapped path. The Louvain equations,
row visitation, summation order, seed, convergence criteria and work allowance
are unchanged. Sequential offset and scatter-bucket readers also keep their
16 MiB windows. This repairs measured mapping overhead while preserving the
bounded edge working set.

## Evidence behind the change

The original full HIRISA clustering process was alive when sampled. In the
one-second sample, 505 of 512 active-thread stacks were in mapping/unmapping
calls under the graph row reader. Its original executable and run were preserved.

A standalone probe compiled the exact old mapped-reader definition and the new
buffered reader. Both read all **34,707,084** original HIRISA edge records and
agreed on every row, column and FP64 bit. The input edge/offset hashes match the
already-qualified graph receipt. A fixed seed-7 shuffle then selected 131,072
rows containing 2,821,220 edges, with mapped/buffered/buffered/mapped trial order.

| Read-only primitive | Mapped seconds | Buffered seconds |
| --- | ---: | ---: |
| Published graph inode | 1.022–1.042 | 0.106–0.108 |
| Active run's private snapshot inode | 1.197–1.516 | 0.126–0.131 |

Each buffered trial filled 141,559 pages and read 579,825,664 bytes, including
rereads. These are bytes delivered to the buffer, not physical-disk traffic.
The old graph calculation and other qualification work ran concurrently. These
observations guide the storage change; they are not a controlled end-to-end
speed comparison or a completed optimized million-cell clustering result.

## Native regression and lifecycle qualification

The final release executable is
`0213c1ca971f1acb0128f6c330d2f1aea10e47e452912a8dd965e38f3ccfda15`.
Its environment binds all 533 authored production source files, compiler, HDF5,
hardware and actual macOS identity. Twelve tests in three suites and the complete
release build pass. The tests cover exact records across pages, a file exceeding
one map window, malformed tails, changed files, symlinks, axis/work rejection,
stable aggregation and exact resident/file solver agreement.

Twenty-two graph fixture commands and 22 clustering fixture commands pass,
including six expected rejections in each suite. Fitted/query, exact/HNSW,
replay, the frozen solver oracle and rehashed corruption controls are covered.
Fresh complete Baron (8,569 cells) and Hagai (13,863 cells) graph records match
the earlier qualified bytes. Their entire clustering result JSON is byte exact
against the frozen results, and native reconstruction passes. Both small edge
files remain in one mapping, so their buffered-read counters are correctly zero.

The execution report adds optional `edgeBufferLoads`, `edgeBytesRead` and
`maximumEdgeBufferBytes`. Load/byte accounting describes deterministic complete
page fills, independently of interrupted or short OS reads. It excludes mapped
small graphs, offsets, scatter buckets and parent reconstruction. Existing
edge-visit accounting continues to cover algorithm work.

The first always-buffered candidate, its source snapshot and successful tests
are retained separately. The final selection keeps the efficient single-map
path for small graphs and aggregate levels. No old receipt was rehashed to claim
a new execution. This archive does not contain a completed million-cell partition.

## Separate integration storage prototype

The unchanged production integration matrix was exercised at 1,612,594 rows and
100 columns with a synthetic storage pattern, not biological observations. Its
private scratch file occupied 1,651,296,256 bytes with one 64 MiB map. The original
5% block contains 80,629 shuffled rows. Sorted gather/scatter restores original
row order before arithmetic; all selected values and before/after arithmetic bits
agree across the following paths.

| Storage-only read/update pattern | Seconds in repeated trials | Buffered value bytes |
| --- | ---: | ---: |
| Original shuffled access | See retained per-trial logs, approximately 3.4–4.3 | No row batch |
| Gather complete block | Approximately 0.13–0.14 | 64,503,200 |
| Fixed 8,192-row tiles | 0.358–0.372 | 6,553,600 per tile |

The fixed tile avoids growing a buffer with the full block. Removal, update and
addition remain separate phases, each in the original logical row order. Scratch
was removed after each probe. This is a design candidate for the next integration
storage change; the production solver does not yet use it. Full solver/oracle
comparison, complete HIRISA execution and biological-preservation gates remain
required. Existing negative preservation evidence remains unchanged.

## Retained archive

The [archive manifest](evidence/2026-09-10-storage-access/manifest.json) covers
1,509 members, including both clustering candidates, source snapshots, probes,
all small-cohort outputs, fixtures, controls and runtime identities. Compressed
logical payload is 85,985,455 bytes; identical chunks reuse earlier graph data.
Manifest SHA256:
`ddecc4f4afcd469fc517e919f6693b1443b4d4b0d255d4851a840bd297ea4694`.
Verify every encoded member and full decoded source with
`python verify_archive.py evidence/2026-09-10-storage-access`. Large original
H5AD files and executables remain external with their identities recorded.
