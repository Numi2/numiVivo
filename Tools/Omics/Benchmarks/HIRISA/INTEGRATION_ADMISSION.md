# Complete-cohort integration witness admission

The file-backed ridge solver admits 1,612,594 HIRISA cells with 100 clusters,
20 PCs and ten iterations under its explicit work allowance. Its membership
output needs 2,580,150,400 bytes. Publication could write that output, but the
independent 1,600,000,000-byte snapshot/replay ceiling would reject it later.

Shared witness bounds now follow the supported two-million-cell axes and the
unchanged method maxima of 100 clusters or neighbors. Membership and MNN anchor
limits are both 3,200,000,000 bytes. MNN retains at most k upper-level candidates
per row; mutual anchors are a subset of those candidates across all levels.
Publication, snapshot and replay use the same anchor bound. Method defaults,
solver equations, row order, work limits and resident-memory allowances are unchanged.

## Verification

Thirteen native tests in four suites pass on the physical M4 Pro Mac mini.
A sparse 1,600,000,016-byte witness exercises real snapshot copying and hashing
above the old ceiling; the clone fingerprint agrees, and a file exceeding the
new bound is rejected. Full file/resident solver trajectories agree across three
seeds and mapping windows. The complete release build also passes.

The frozen release executable is
`a4ad3dae32dc3568b4b9832d2376a5888af49249881abc907fbb40921f4973a3`.
Its environment manifest records all 532 authored production source hashes,
compiler, HDF5 and hardware. Fresh self-contained CLI fixtures run on the laptop
with that host's actual OS identity. Existing Mac mini receipts are not relabeled.

- 60 ridge lifecycle commands pass, including 16 expected rejections. All frozen
  legacy matrix bits agree; fitted/query, donor/batch, fixed/adaptive modes,
  exact/HNSW downstream graph, clustering and embedding are covered.
- 45 MNN lifecycle commands pass, including 11 expected rejections. Scanorama
  1.7.4 independently reproduces fitted/query anchors exactly and corrected
  coordinates within 9.658940314238862e-15 maximum absolute error.
- Rehashed corrupted witnesses/parents, unsupported methods and invalid resource
  plans are rejected; scratch cleanup and exact repeat receipts pass.

The [evidence archive](evidence/2026-09-10-integration-admission/manifest.json)
contains 1,353 members (1,538,188 compressed bytes), including full fixture
outputs, controls, source/runtime identities, logs and the original admission
audit. Manifest SHA256:
`c1b8c9c5bd356b3551bc8d642c0af4e850b15e004e3ef31737e96d505231a78d`.
Verify with `python verify_archive.py evidence/2026-09-10-integration-admission`.

## Remaining full-cohort work

These checks repair artifact compatibility; they do not qualify complete HIRISA
integration or biological preservation. Full source PCA must be refitted under
the new executable before its native integration receipts can be used.

The exact MNN donor-matching workload alone requires 20,684,720,015,420 scalar
terms, beyond the supported ten-trillion ceiling before anchor assembly. Its
k=20 base latent estimate is 3,225,188,000 bytes before anchor buffers and source
metadata. Raising the artifact ceiling does not solve those algorithmic limits.

The ridge baseline also needs a storage-access measurement before full execution:
shuffled sweeps traverse a 1,651,296,256-byte membership scratch file through one
64 MiB mapping. Frequent remapping is a source-derived concern, not a measured
performance result. Preserve numerical order when qualifying storage changes.
Earlier negative cell-type preservation results remain valid; this repair does
not turn them into biological success.
