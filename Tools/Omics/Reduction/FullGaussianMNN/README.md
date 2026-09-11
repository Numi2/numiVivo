# Full Gaussian correction recovers signal, approximate matching still fails

Keeping the earlier eligible-level HNSW matcher fixed and restoring full-anchor
Gaussian smoothing recovers the measured Ding preservation gates and most Kang
gates. **Kang megakaryocyte recall still fails its original margin.** This fixed
candidate remains experimental; the native production API keeps exact matching.

Every original Hagai, Kang and Ding cell participates. This is development on
inspected cohorts, not independent biological validation or prospective prediction.
No labels or programs enter fitting, no cohort is removed, and no parameter is
tuned after the result. See the predeclared [protocol](PROTOCOL.md).

## Complete original cohorts

| Cohort | Cells | Exact-anchor precision / recall | Fit / replay seconds | Original biological gates |
| --- | ---: | ---: | ---: | --- |
| Hagai | 13,863 | 99.896% / 99.862% | 2.843 / 2.839 | All applicable measured gates pass |
| Kang | 24,673 | 99.994% / 99.948% | 13.780 / 13.798 | Megakaryocyte recall fails; four classifier strata remain unavailable |
| Ding | 44,031 | 99.748% / 99.014% | 28.513 / 28.509 | All measured gates pass; original labels remain partial |

Kang overall type balanced accuracy improves from 0.876995 in the earlier local64
trial to 0.898019, meeting the fixed allowance against the uncorrected baseline
0.910778. Within-stratum program Spearman improves from 0.478409 to 0.553630,
also meeting its original margin. But megakaryocyte recall is 0.062218 below the
uncorrected baseline, exceeding the fixed 0.05 allowance. The exact production
MNN meets that same gate. The four missing condition-classifier strata remain
missing; they are not counted as passing tests.

Ding type balanced accuracy improves from 0.624435 to 0.658914, compared with
0.660902 before correction. All per-type recall and global/within-stratum program
margins now pass. Only 29,411 of the original 44,031 cells carry source labels;
all cells remain in fitting and every available original label is evaluated.

The matching neighbors and mutual anchor arrays are identical to the previous
[eligible-level local64 trial](../LocalMNN/README.md). Thus this comparison
isolates the smoothing change while holding the approximate matcher fixed.
It shows that local smoothing explains much of the previous measured loss, but
that retaining over 99% of exact anchors still does not guarantee all biological
preservation gates. This is a result for these fixed development settings, not
a proof that every approximate matching method must fail.

## Numerical method and execution scope

The driver reuses the exact previously qualified native HNSW library and the
production [tiled Gaussian library](../GaussianKernel/README.md), verified by
source and executable hashes. HNSW remains M=16, construction ef=200, query ef=128,
seed=7, with separate eligible-prefix/suffix indexes and FP64 Manhattan matching.
All mutual anchors and duplicate records remain in the full Gaussian sum.

The research assembly reuses the existing Python panorama owner through an
optional kernel callback. Its default local64 route reproduces every original
Hagai witness and report array exactly after this refactor. The new callback uses
32-query/256-anchor native tiles, direct squared distances, sigma=15 and the
qualified scalar underflow fallback. The driver admits twice the exact distance
work per correction and records actual terms, including any fallback work.
No dense cells by genes array is constructed.

Native matching took 1.626, 3.753 and 7.262 seconds respectively. Entire fit/replay
execution took 99.99 seconds, with peak RSS 251,478,016 bytes. Reported fit times
exclude writing/hashing full diagnostic witnesses; the process total includes
that work. The physical M4 Pro's original mapped-reader clustering baseline
remained active. These research-driver timings do not establish end-to-end
speedup over the production Swift API, scverse or million-cell integration.

## Evidence and reproduction

All complete score, anchor, neighbor and per-step snapshot arrays replay exactly.
The biological evaluators use the unchanged original input NPZs, thirty-neighbor
protocol, folds, classifiers, programs and margins. Full neighbor/prediction
arrays and all failures are retained with output identities frozen before scoring.

Use the original cohort manifest and evaluation spec, adapting paths without
changing source identities:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python \
  Tools/Omics/Reduction/FullGaussianMNN/run.py \
  --source ORIGINAL_COHORT_ROOT --manifest ORIGINAL_MANIFEST.json \
  --library QUALIFIED_LOCAL_NEIGHBORS.dylib \
  --gaussian-library QUALIFIED_GAUSSIAN.dylib \
  --protocol Tools/Omics/Reduction/FullGaussianMNN/PROTOCOL.md --out NEW_FIT
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python \
  Tools/Omics/Reduction/FullGaussianMNN/check_numerics.py \
  --source ORIGINAL_COHORT_ROOT --root NEW_FIT \
  --previous PRIOR_ELIGIBLE_LEVEL_FIT --out NUMERICAL.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python \
  Tools/Omics/Reduction/FullGaussianMNN/evaluate.py \
  --source ORIGINAL_COHORT_ROOT --root NEW_FIT \
  --spec ORIGINAL_EVALUATION_SPEC.json --out NEW_EVALUATION
```

Build the neighbor library with the recorded [LocalMNN command](../LocalMNN/README.md)
and the Gaussian shared library from `Sources/NumiVivoCore/OmicsGaussian.cpp` with
C++23, `-O3 -shared -fPIC`, the Core include path and `-framework Accelerate`.
Environment receipts bind the actual reused libraries and all vendor sources.
The numerical checker independently evaluates every full Gaussian update with
NumPy direct coordinate differences, exponentials and weighted sums. A completed
fit is never restarted just because an observation expires.

Independent checks completed all 650,248 correction-row updates. Maximum
Gaussian delta discrepancy was 8.42e-14, all recorded matching distances
were exact, and final coordinates reconstructed exactly. All alignment orders
match the original exact MNN; relative coordinate Frobenius differences are
0.007865, 0.024625, 0.127885 (Hagai, Kang, Ding). Numerical agreement and high
anchor recall do not override the failed Kang biological gate.

The [verified archive](../evidence/2026-09-11-full-gaussian-mnn) contains 49
members and 51,558,445 stored bytes; manifest SHA-256
`4b1cfe9f658268aaf8a4cc58ee86611771a39842c1ab6ce95b5c2c2ea0a873a5`.
It includes complete corrected scores, matching/mutual anchors, biological
neighbor/prediction arrays, reports, numerical checks and source. All 34 full
step-witness files remain externally retained with exact paths and hashes.
After archive verification, identical remote witnesses were removed only after
local-copy and open-handle checks. The [storage receipt](../evidence/2026-09-11-full-gaussian-mnn-storage.json)
records current laptop placement and 302,608,384 recovered bytes. The archive
predates this relocation; restore these witnesses from the receipt's retained
root before replaying checks on the Mac mini. Original cohorts, binaries,
score/anchor results and active clustering baseline remain intact.

This fixed hypothesis is complete and should not be rerun or tuned to erase its
failure. Further scaling work must preserve biological gates, including rare
labels; neither near-exact anchor recall nor recovered overall accuracy is enough.
