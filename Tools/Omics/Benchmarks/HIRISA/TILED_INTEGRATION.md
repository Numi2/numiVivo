# Bounded integration row batches — 10 September 2026

The production integration solver now gathers and scatters shuffled rows in
batches of at most 8,192 when its latent matrices exceed one mapping window.
It retains the original logical arithmetic order and all three whole-block
phases: remove old memberships, compute replacements, add new memberships.
Resident and single-window matrices retain direct row access. Defaults, ridge
models, stopping rules and work admission are unchanged.

At 100 clusters the tile value payload is at most 6,553,600 bytes (6.25 MiB),
independent of cohort size. Array/index overhead, transient per-row probability
vectors and mapped pages are additional memory. The report's optional
`maximumBufferedValueBytes` describes this tile payload bound, not process RSS.
Whole-batch validation rejects oversized batches, duplicate or out-of-range
indices, mismatched widths and nonfinite values before writing. Private scratch
may be interrupted by cancellation or IO failure; final publication is transactional.

## Verification

The complete native release build passed. All 533 production source hashes are
bound to executable SHA256
`c819027a4f5f723ae2c5f4c04edd21dea2888f33b925c2d04763d00d56fc1c47`.
Native execution used macOS 26.6, build 25G72. Fresh laptop fixtures used the
same executable with the laptop's actual build 25G5028f; foreign receipts were
not relabeled.

- 14 tests in four suites passed. These include exact file/resident trajectories
  for three seeds, signed-zero/subnormal batch values, full/tail batches,
  invalid-batch rejection without writes, and the large witness admission check.
- A separate scoped solver check uses 163,881 synthetic rows, three clusters and
  a 16 KiB test mapping. Each normal block contains 8,194 rows and crosses the
  fixed 8,192-row tile boundary. Current file execution matches the exact prior
  resident solver in every matrix record and all scalar diagnostic bits, for
  fixed ridge with seed 7 and expected-mass ridge with seed 19. Production solver
  and matrix files are compiled unchanged; scoped dependency support and the
  unused PCA overload's type adapter are retained explicitly.
- All six complete Kang/Hagai trajectories (24,673/13,863 cells; seeds 7, 19, 41)
  reproduce the previously published frozen-oracle-qualified matrices and
  numerical diagnostics. Fresh PCA ancestors match every checked PCA payload.
  All six native replays pass; 14 native commands in total. These smaller
  matrices fit one window and exercise the direct-access path.
- 60 fitted/query ridge lifecycle commands pass with 16 expected rejections and
  exact frozen matrix bits. Another 45 MNN commands pass with 11 expected
  rejections; Scanorama anchors match and maximum coordinate error is
  9.658940314238862e-15. Downstream exact/HNSW graphs, clustering, embeddings,
  repeated receipts and rehashed corruption controls are included.

## Full-shape storage measurement

The unchanged production matrix and its new batch APIs were measured with
1,612,594 synthetic rows, 100 columns, 1,651,296,256 scratch bytes, one 64 MiB
mapping and 80,629 selected shuffled rows. Every selected value and both
order-dependent sum bit patterns match across all trials.

| Access | Trial seconds |
|---|---:|
| Shuffled | 3.4805 / 4.1347 |
| Batches of at most 8,192 | 0.3195 / 0.3114 |

The whole storage driver used 73,007,104 maximum resident bytes and removed its
scratch. These are storage-only measurements under the recorded concurrent
workload, not complete integration or a controlled end-to-end speed comparison.

## Remaining full-cohort work

Complete HIRISA integration and biological preservation are still unqualified.
The prior Kang NK-cell recall failure and unavailable evaluation strata remain
unchanged. This storage repair does not improve those biological outcomes.

The subsequent [scratch lifetime repair](SCRATCH_LIFETIME.md) releases
solver-only scratch after solving and each output matrix after successful
serialization, with cleanup/replay qualification. Full execution then needs fresh
runtime-bound PCA ancestors and the original cohort, followed by independent
numerical and biological-preservation checks. Exact MNN scaling and Metal work
remain separate.

## Archive

[Manifest](evidence/2026-09-10-integration-tiled/manifest.json): 1,539 members,
263,167,591 logical compressed bytes. Identical prior compressed payloads are
reused, including 64 members from the earlier full integration archive, so these
are not all new Git objects. Large original H5AD inputs and executables remain
external with verified paths and identities.

Manifest SHA256:
`02e97fd1bc7f8d0283d093a04d4d314e09d736c962f9c3c40c6a5b611141d9b0`.

Verify with
`python verify_archive.py evidence/2026-09-10-integration-tiled`.
