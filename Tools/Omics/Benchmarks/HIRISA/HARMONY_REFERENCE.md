# Complete HIRISA independent integration reference

A subsequent [within-library decoder development experiment](PROGRAM_CALIBRATION.md)
meets the unchanged program-loss margins with 28/32 sensitive controls, four
insufficient. It identifies objective dependence without replacing the original
failures or establishing independent biological preservation.

The subsequent [full-cohort program diagnostic](INTEGRATION_PROGRAMS.md) fails
three control-sensitive comparisons for native integration and all three Harmony
seeds, with 18/32 comparisons insufficiently sensitive. Earlier coarse results
remain valid at their stated scope; complete biological preservation is unqualified.

All three Harmony 2.0.0 reference fits and response evaluations are complete for
**1,612,594 cells** using
the same original 20-PC matrix and donor metadata as native integration. No
cells, genes or metadata levels are selected from the native response results.
The reference protocol was frozen after inspecting the native seed-7 response
diagnostics; this timing is explicit and is not a preregistration claim.

The settings are the unchanged prior Kang/Hagai reference settings: seeds
7/19/41, 100 clusters, theta 2, fixed lambda 1, sigma 0.1, tau 0, block fraction
0.05, ten Harmony iterations, four k-means iterations, cluster tolerance 0.001,
Harmony tolerance 0.01, alpha 0.2, batch proportion cutoff 1e-5 and one BLAS
thread. Independent initialization differs from native initialization.

`run_harmony_integration.py` verifies every input record and metadata assignment,
loads only the latent matrix into the resident reference, and writes corrected
coordinates in the same explicit row/column/Float64 record format used by the
frozen response evaluator. It independently reads every written coordinate back
and checks exact bits against the reference output. It records full source,
protocol, package code/binary and output identities, objective history, timing
and peak process RSS. Harmony may allocate a resident cells × clusters
workspace; this is a reference measurement, not a claim of out-of-core Harmony
or a resident native integration implementation.

The wrapper is qualified on **all 24,673 Kang cells**: seed-7 output agrees
exactly, with maximum absolute discrepancy zero, with the previously retained
complete Harmony 2.0.0 result. This checks metadata coding, matrix orientation,
the pinned options and output transport on experimental data.

The frozen coordinator waits for the existing native full replay to pass and
at least 8 GiB free disk before starting. Each seed then requires at least 4 GiB
free before fitting, runs sequentially, and uses the same frozen 79-fold
response evaluator as native integration. Failures and all four insensitive
classification controls remain visible. The existing native seed-7 run is not
restarted or replaced. Completion/results are recorded separately from this
execution specification.

## Complete comparison

Native seed 7 and all three independent reference seeds pass the frozen coarse
response-preservation margins. The exact baseline and all 79 matched folds are
identical across methods; classification sensitivity remains insufficient in
the same four contrasts.

| Method/seed | Largest contrast-mean accuracy loss | Largest fold loss | Largest contrast-mean response drift |
| --- | ---: | ---: | ---: |
| Native / 7 | 0.002545 | 0.019370 | 0.138859 |
| Harmony / 7 | 0.002527 | 0.020011 | 0.137067 |
| Harmony / 19 | 0.002387 | 0.020353 | 0.135625 |
| Harmony / 41 | 0.002384 | 0.019078 | 0.134844 |

Reference fit observations are 85.65/63.96/96.12 seconds for seeds 7/19/41,
with maximum process RSS 3,861,790,720 bytes. These fit-only times exclude native
parent replay/PCA and are not directly comparable to the native end-to-end
publication/replay timings. The experiments ran on a concurrently used Mac mini.
No speed advantage is inferred from those unlike timing scopes.

The [reference report archive](evidence/2026-09-10-harmony-reference) contains
39 members/204,065 compressed bytes, manifest SHA256
`907050e45ce611df147a3bad644917328f70f90268b0b45c67fddf6d212919c1`.
It retains the frozen protocol, coordinator/runner, complete Kang wrapper
qualification, all fit reports, package identities, 79-fold response results
for every seed and terminal logs. All three full corrected matrices remain
preserved externally with exact rechecked hashes. `archive_harmony_reference.py`
validates completion and score/report bindings before archiving.

Full-cohort preprocessing is transductive. Three reference seeds versus one
native seed do not qualify native initialization robustness. A method comparison
on the coarse response metrics does not establish marker/program gradients,
rare-cell preservation, learned reference annotation, prospective prediction,
causal validation or a controlled performance advantage.
