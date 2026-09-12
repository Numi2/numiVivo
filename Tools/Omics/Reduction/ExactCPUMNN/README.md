# Exhaustive native CPU matching preserves the exact benchmark

The C++ exhaustive matcher recovers **every native exact anchor** across all
**82,567 original Hagai, Kang and Ding cells**, and final corrected coordinates
agree within 1.30e-13. Every evaluable original biological gate passes, including
Kang's condition, per-type recall and within-stratum program gates. Four original
Kang treatment-classifier strata remain unavailable; Ding labels remain partial.
These limitations are retained rather than counted as passes.

| Cohort | Exact anchors | Matching distance evaluations | Matching seconds | Fit seconds |
| --- | ---: | ---: | ---: | ---: |
| Hagai | 49,893 | 126,392,182 | 0.631 | 1.809 |
| Kang | 159,513 | 515,933,916 | 2.799 | 12.489 |
| Ding | 131,282 | 1,600,487,330 | 6.933 | 26.859 |

This reuses the native exhaustive kernel qualified in the distance-tail trial,
but executes it for every eligible query; HNSW and tail selection are not called.
Original normalized PCA, Manhattan distance, k=20, grouped-source-index ties,
full-Gaussian sigma=15 and panorama assembly stay unchanged. Work is partitioned
by donor/batch level and search direction, with a 500-million distance budget
per partition and all partitions required. The measured total is independently
checked against the complete cross-level directed pair count. There is no
cells-by-genes dense matrix, approximate matching or cell exclusion.

## Evidence and scope

Actual fitting and full replay pass for all cohorts. Every selected distance,
all mutual anchor identities and every final coordinate are compared with the
original native exact result. Assembly order matches. Fresh biological scoring
uses the unchanged original cohorts, baselines, folds, programs and margins.
The earlier wider-query and distance-tail failures remain available; they are
not repaired by changing thresholds.

The native library is reused with verified binary, source and header identities;
only the Python research driver's matching invocation changes. Fit-time source,
protocol and inputs were frozen before execution. This is an implementation and
development benchmark, not a new independent biological validation. The product
Swift owner has not yet been switched to this kernel. Million-cell MNN execution,
end-to-end product memory/timing and independent preservation remain unqualified.

Observed fit times are lower than the earlier width-512/hybrid runs, but these
are separate sequential runs, not a randomized paired performance study. They
exclude witness serialization/hashing and do not establish application or scverse
speedup. They justify testing a production integration of the exact kernel
before accepting approximate matching's measured biological trade-offs.

## Reproduction

`python3 verify.py` checks the archive, source/library identities, output/array
bindings, complete cohort accounting and all evaluable gates. It explicitly
retains the missing Kang classifier flag. [summary.json](summary.json) records
all gates and timings; [manifest.json](manifest.json) binds every full NPZ
matching, Gaussian and biological witness retained on the Mac mini. The archive
includes code, headers/library, frozen protocol, commands, reports and logs.
Original cohort artifacts, evaluator inputs and the qualified Gaussian library
are external dependencies; the result archive is not a standalone raw-data replay.

Restore those exact dependencies and run archived `FullGaussianMNN/run.py` using
`PROTOCOL_EXACT.md`, the recorded native/Gaussian libraries and a fresh output
path. Run `check_exact.py` and `LocalMNN/evaluate.py` with the retained evaluation
spec afterward. Adjust absolute paths without changing identities. The exact
kernel remains quadratic; partitioning is a work-accounting boundary, not an
asymptotic scaling improvement.
