# Metal kNN distance blocks: real-Hagai research benchmark

A Metal FP32 squared-Euclidean distance kernel reduces median process time by
**35.1% (1.54x)** versus the same harness's scalar Swift FP32 CPU path, after
moving the Metal output-buffer pointer lookup outside the selection loop.
All **13,863 Hagai cells**, 20 original PCs and 20 neighbors including self are
retained. All 192,168,906 directed non-self candidates are evaluated. Neighbor
selection remains on CPU; every candidate is searched, using FP32 arithmetic.

| Version | CPU elapsed seconds | Metal elapsed seconds | Median CPU / Metal |
| --- | --- | --- | --- |
| Initial repeated pointer lookup | 3.026, 2.650, 2.648 | 3.536, 3.513, 3.380 | 2.650 / 3.513 |
| Hoisted pointer lookup | 2.974, 2.572, 2.589 | 1.708, 1.682, 1.674 | 2.589 / 1.682 |

Both fixed CPU/Metal/Metal/CPU/CPU/Metal sequences are retained. The initial
slower GPU candidate was not discarded. These are separate process timings on
the physical Apple M4 Pro: loading/validating score records, FP32 conversion,
shader compilation, distance computation, CPU neighbor selection and file output
are included. Source SHA checks in the outer driver are outside these timings.
Caches are exercised; no cold-cache, randomized confidence interval, Scanpy
speedup or full product pipeline comparison is claimed. The final CPU peak RSS
is about 20.7 MB, Metal 28.7–29.1 MB. GPU-command time is about 0.136 s; it is not
the complete elapsed time.

## Numerical results and admission boundary

The independent SciPy FP64 checker searches the complete candidate axis with
source-row tie breaking. Both CPU and Metal retain **100% neighbor membership**.
One row changes neighbor order versus FP64 in both paths. CPU and Metal indices
are identical, but their distances are not bitwise identical. Maximum scaled
squared-distance error is 5.83e-7 CPU and 6.04e-7 Metal (scale=max(1, reference)).
All three repetitions per backend produce identical indices and distances.
The same numerical findings hold for the initial candidate.

FP32 precision and graph order matter downstream. This result does not establish
unchanged graph weights, embeddings, clustering, batch mixing or biological
preservation. It is a research harness, **not an installed native kNN backend**.
It does not replace the native FP64 default. It retains the full resident PCA
score array and a 64-query-by-N distance tile (3,548,928 bytes for this cohort),
not an N-by-N matrix or a dense cells-by-genes array. Million-cell out-of-core
execution is not qualified. The harness rejects sizes outside 20–50,000 cells
and 1–64 components; broader admission/cancellation/rejection coverage is still
required before owner integration.

The next gate is integration into the native tiled neighbor owner with explicit
FP32 provenance, cancellation and bounded buffers, followed by complete-cohort
CPU/scverse comparisons and downstream biological-preservation checks. A single
cohort's timing gain cannot justify changing the default.

## Reproduction and evidence

Compile `Main.swift` with `swiftc -O -swift-version 6 -framework Metal` on a
physical Apple GPU. `run.py` records source/binary/script identities before the
six runs; its recorded source paths require the retained Hagai PCA bundle.
`check.py` verifies all score coordinates, repetition equality, complete FP64
neighbor selection and every selected distance. No cells are subsampled.
`verify.py` checks archived evidence without rerunning the experiment.
The archive retains initial and final sources, protocols, logs and numerical
reports. Executables and all output arrays remain externally hash-bound in its
manifest. The compile deprecation warning is retained.
