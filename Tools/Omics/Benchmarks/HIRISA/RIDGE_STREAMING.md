# Sequential ridge passes at complete HIRISA shape

File-backed integration now accumulates all active cluster/level ridge sums in
one ascending-cell pass and applies the fitted effects in a second ascending-cell
pass. The previous implementation scanned the full latent and membership
matrices separately for each cluster, and revisited score rows for each active
cluster/level correction. This change removes repeated storage reads and row
allocations without changing the estimator or arithmetic order within a result.

For each cluster/level/component, the accumulator still receives exactly the
same `weight * score` products in ascending cell order. Fits run in ascending
cluster order using the unchanged fixed/adaptive Schur solver and active-level
threshold. Each cell then receives subtractions in exactly ascending cluster
order. Inactive cluster/level effects are skipped rather than represented as
zero-valued operations. Original PCA values remain the correction target.
Initialization, shuffled block updates, objectives and stopping rules are unchanged.

The new branch applies when any relevant file matrix spans multiple mapping
windows. Resident and single-window cases retain the prior direct path.
The weighted-sum table contains at most 100 clusters × 128 levels × 64 components
× 8 bytes = **6,553,600 value bytes**. The effect table has the same maximum
value payload. These tables are independent of cell count; arrays, optionals,
indices, per-row temporaries, mapped pages and existing metadata add memory.
At HIRISA's 100 clusters, five donors and 20 PCs, each table has only 80,000
value bytes. There is no full cells × clusters or cells × genes resident cache.
Both passes check cancellation every 2,048 rows; malformed axes/effect tables
are rejected before score writes. Publication retains ownership of scratch
cleanup on failure or cancellation.

## Verification

The release executable SHA256 is
`67343146e05b783b1bb0806cb765a536890109a6b50816d34617144e0af7c114`.
All 533 production source identities are frozen in the archived environment.
Nineteen tests in four suites pass, and the complete release build passes.
The additional checks exercise exact scalar order across mapping boundaries,
inactive effects, signed zero, validation before writes and cancellation.

All six full Kang/Hagai seed-7/19/41 native outputs retain exact PCA, membership,
assignment and corrected-score bytes, complete trajectory diagnostics, and
successful native replay. Those cohorts fit the usual 64 MiB window and test
the direct path. A separate scoped build uses the unchanged production solver
with a **16 KiB scratch window** to force the new path on all six complete
cohorts. Every matrix and scalar diagnostic bit matches the original frozen
native results. This scoped check supplements the native receipts; it does not
rewrite their executable or host bindings.

The first scoped harness rejected the metadata's `group` field when decoding it
as a strict bare identity, before fitting. Its failure is retained. The corrected
harness explicitly extracts sample ID/barcode, matching native publication;
production code, cohort, options and reference outputs did not change.

A separate 163,881-cell synthetic trajectory compares the new file solver with
the exact preceding resident solver for fixed seed 7 and adaptive seed 19. Each
8,194-row shuffled block crosses the 8,192-row tile boundary. Every output and
diagnostic bit remains exact. The 60 fitted/query ridge commands (16 expected
rejections) and 45 MNN commands (11 expected rejections) also pass; Scanorama
anchors and the prior coordinate tolerance remain satisfied.

## Full-shape phase measurement

The local probe uses **1,612,594 rows, 20 PCs, 100 clusters and five levels**,
with file-backed inputs and two separate output matrices. It checks every
weighted statistic and every corrected scalar bit against the previous loop
order and removes all generated scratch after success.

| Phase | Previous passes | Sequential passes |
| --- | ---: | ---: |
| Weighted statistics | 34.1667 s | 2.0413 s |
| Apply effects | 91.9722 s | 1.9081 s |

These are sequential observations from a synthetic phase probe on a concurrently
used MacBook. They exclude initialization, membership optimization, PCA,
publication and replay. They are not controlled end-to-end or CPU/Metal speed
comparisons. The full HIRISA publication using the earlier frozen executable
is a separate experiment; this new executable has not yet completed a full
HIRISA publication/replay.

The [evidence archive](evidence/2026-09-10-integration-ridge-stream) contains
1,587 members and 430,668,961 logical compressed bytes, including both scoped
harness attempts, exact previous source, complete cohort outputs and checks,
fixtures, tests, logs and phase measurements. Existing verified compressed
payloads supplied 430,202,053 bytes; only 466,908 new compressed bytes were
transferred. Its manifest SHA256 is
`a0b00b440ee6759b8e5873f65b03c07dcaaa29f8475657d4105e47fca4788737`.
`archive_ridge_integration.py` verifies prerequisite outcomes and can reuse
existing archives without copying unchanged payloads. `verify_archive.py`
checks every stored/decoded member and full source reconstruction.

This storage change does not repair known Kang NK-cell preservation failures,
qualify rare-cell or marker/program preservation, or complete the broader goal.
