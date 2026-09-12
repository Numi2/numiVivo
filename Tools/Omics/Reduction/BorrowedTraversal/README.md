# Bounded sparse PCA arithmetic

The native sparse PCA owner now borrows its input and result buffers once per
project/transpose pass and inlines record traversal. The reader checks its mapping
once per 16 MiB window. No dense cells × genes matrix is introduced. Private cache
coordinates and values remain validated once before fitting. Source order and FP64
multiply/add order are retained for CSR and CSC input; cancellation is checked at
window boundaries. Existing random record access is unchanged.

This follows the [unpromoted mapping-only experiment](../WindowTraversal/README.md).
The combined borrowed-buffer and inlining change is measured below; these runs do
not isolate the contribution of each compiler optimization.

## Full real-data execution

Physical M4 Pro, 2026-09-12. Actual native H5AD PCA publication, source reconstruction
and output writing included, subsequent hashes excluded. Three runs each in fixed
old/new/new/old/old/new order per cohort, one-worker environment. Baseline binary
is the qualified native Metal kNN CLI; PCA source was unchanged through `695cf43a`.

| Cohort | Cells | Original median | New median | Lower elapsed |
| --- | ---: | ---: | ---: | ---: |
| Kang | 24,673 | 4.068 s | 3.742 s | 8.0% |
| Hagai | 13,863 | 7.386 s | 6.568 s | 11.1% |

All twelve commands succeeded and every score, loading and metadata artifact
matched the qualified cohort output byte-for-byte. The first new Kang run took
4.095 s, slower than the old median; it remains included. This small warm-host
comparison has no confidence interval or cross-device performance claim.

Separate single Kang `/usr/bin/time -l` runs measured peak RSS of 139,083,776 bytes
original and 167,985,152 bytes new, with elapsed 4.02 and 3.69 seconds. The new run
had higher observed peak memory; no memory reduction is claimed. These are separate
observations, not the six-run timing medians. The mapping bound remains 16 MiB;
metadata, results and solver basis have their existing resident-memory limits.

## Validation and reproduction

`Controls.swift` checks two full traversals after prior random access, a partial
final window, unsorted coordinate project/transpose against sequential FP64 sums,
visit accounting, and cancellation before traversal. These pass against the actual
scoped library. AddressSanitizer also passes the extracted, unchanged reader and
reduction cache owners with only error and count-limit stubs; this is focused
memory-safety evidence, not an instrumented full application or full package test.

The scoped CLI build and controls use `validate.sh`; `check.py` and `hagai.py` retain
exact source/plan/output paths. Build with `Tools/Omics/H5AD/build.sh BUILD --with-cli`.
The sanitizer invocation is `swiftc -swift-version 6 -O -sanitize=address
-parse-as-library SanitizedOwners.swift SanitizedControls.swift -o sanitized-controls`,
then `./sanitized-controls asan-controls.bin` in a fresh directory. Numerical benchmark
and control outputs, compiler logs and sanitizer sources are retained here.
`runtime-manifest.json` binds source, plans, executables and all twelve numerical
outputs. Large runtime artifacts remain at `/Users/n/numivivo-pca-borrowed-20260912`.
Run `python3 Tools/Omics/Reduction/BorrowedTraversal/verify.py` to check current owner
hashes and retained numerical/timing records. It does not rerun external datasets.

This is a CPU sparse PCA improvement. GPU PCA, million-cell PCA, a matched scverse
PCA benchmark, and general biological prediction remain unqualified.
