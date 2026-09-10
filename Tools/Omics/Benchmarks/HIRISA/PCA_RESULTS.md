# Complete HIRISA PCA qualification

Native PCA publication and independent comparison of all 32,251,880 scores now
pass on the complete 1,612,594-cell, 18,082-gene source. Native replay also passed
and reconstructed the same receipt from the archived source.
The [execution](PCA_EXECUTION.md) and [comparison](PCA_COMPARISON.md)
specifications separate resource admission from numerical and biological claims.

## Completed measurements and checks

The complete-source measurement selects 2,000 genes and 141,471,392 nonzeros.
The unchanged native COO representation requires 2,263,542,272 cache bytes and
44,704,959,872 entry visits at the frozen 20-component, 128-basis settings.
Only explicit resource allowances change; normalization, selection, arithmetic
and residual gates remain unchanged.

Independent chunked feature moments and pinned Scanpy selection were followed
by SciPy ARPACK over a file-backed centered covariance operator. All source
rows and full-library denominators contribute. The reference fit completed in
85.94 seconds (73.03 seconds for source validation/materialization), using 95
covariance applications including residual checks. Maximum relative residual
was 1.339e-13 and loading orthogonality error was 2.220e-15. Observed peak memory
footprint was 1,420,987,152 bytes; maximum RSS was 3,117,449,216 bytes.
Whole-process time was 86.75 seconds. This is a local reference observation, not a
controlled performance comparison against native execution on the Mac mini.
The fit hash is `40b7d5d64bdae9b87cade4fb36399ab8a8b11bf849586988ad2ba0d4e62646de`.

The updated native owner passes 36 tests in 12 suites, including large record
offsets, score rows past one million, and retention of a 345,933,053-byte QC
payload by the downstream snapshot reader. Both complete Baron and Hagai
regression workflows pass publication, replay, exact repeat and failure gates;
all reduction, metadata and QC fields equal the previously qualified results.
CSR/CSC fixture fitted fields and binary scores match exactly, and canonical
rehashed score/model tampering is rejected by native source reconstruction.

The current full-source native run uses immutable release binary
`d989aaa85157edbb295ad5a50cb8566ecf4321ebafeabe53879f9f0cbd200050`.
Its manifest binds 532 authored source/package files; they match both isolated
checkouts. Complete HIRISA publication passed in 575.29 seconds at maximum RSS
4,331,503,616 bytes. Native replay passed in 580.63 seconds at maximum RSS
4,339,122,176 bytes. Both commands returned zero and produced identical receipts.

Independent comparison checked every original row identity, every QC count and
every score coordinate/value. Metadata and QC bytes exactly equal the previously
verified complete native report. All 2,000 selected genes and all 18,082 feature
statistics pass the frozen gates. The minimum loading subspace cosine is
0.9999999999999973; aligned relative score error is 3.3741e-12 and maximum direct
projection error is 3.9080e-14. Independently measured native loading residual is
1.6307e-11, agreeing with the native report and satisfying the frozen 1e-6 gate.
The comparison took 12.02 seconds at maximum RSS 2,524,184,576 bytes. The native
scores SHA256 is `62e6b0644e2c954749b77946ba748e184c33580614f5756b1bc1b0980df86296`.
These timings come from different hosts and are not a controlled speed comparison.

## Retained failures and interruptions

The unchanged previous runtime rejected the complete source at its 2 GB cache
allowance: exit 65 after 248.33 seconds, maximum RSS 1,215,774,720 bytes.
The independent measurement's first attempt rejected a mistyped expected source
hash; the corrected helper derives that identity from the qualified source
receipt. No source bytes or scientific options changed.

The first large-offset test incorrectly relied on the file offset after
Foundation truncation; explicitly seeking to zero fixed the test setup. The
first native fixture comparison incorrectly required identical private COO
cache hashes across CSR/CSC traversal orders. Every fitted field is exact, while
the raw cache sequence differs. A second harness attempt wrote a noncanonical
tampered receipt, causing receipt-format rejection before reconstruction. Both
attempts are retained; the final control preserves canonical receipt encoding
and reaches the intended reconstruction rejection.

A first updated HIRISA run was explicitly interrupted with SIGTERM after 32.11
seconds, before the expensive fit. Inspection of the existing complete native
report identified another publication cap: quality JSON occupies 345,933,053
bytes, exceeding the old 256 MiB allowance. Metadata occupies 149,506,727 bytes.
The repaired runtime shares a 512 MiB QC bound across publication, replay and
PCA downstream snapshots. The interruption and private staging are retained;
it is not recorded as a native algorithm failure or successful publication.

## Outstanding scope

Metadata, QC and fitted score arrays remain resident.
The C++ HNSW owner still has its separate one-million-row admission limit.
Graph recall, clustering, biological preservation under integration, unseen
perturbation prediction and Metal performance need their own qualification.


## Reproducibility and archive

The [complete evidence archive](evidence/2026-09-10-pca/manifest.json) retains all
score records, loadings, metadata, QC, models, plans, reference fits, command
outcomes, runtime identities and failed attempts. Files larger than 64 MiB are
stored as contiguous gzip parts. Verify both each decoded part and the complete
reconstructed source hash with:

```sh
python Tools/Omics/Benchmarks/HIRISA/verify_archive.py \
  Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-10-pca
```

For reconstruction, group manifest records by `sourcePath`, sort by
`sourceOffset`, and concatenate each gzip-decoded member. `fullSources` records
the complete file sizes and hashes. Exact original H5AD files, native executables
and reconstructible reference CSR scratch remain external, with hashes in the
archive. Public runtime manifests omit device identifiers while retaining the
original private manifest hash.

Rerun the pinned full-source reference with `reference_pca.py --root STUDY`, then
compare a freshly published bundle using `compare_pca.py --root STUDY --bundle
BUNDLE --out NEW_RESULT`. Both preserve all source rows and the full-library
normalization denominators. Native publication/replay use the ordinary
`singlecell-h5ad-pca` and `singlecell-h5ad-pca-verify` product commands with the
archived plan and immutable runtime. The publication and comparison scripts
refuse to silently overwrite earlier outcomes.
