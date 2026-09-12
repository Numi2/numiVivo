# Hagai clustering after Metal FP32 neighbor construction

The actual product clustering path produces **identical CPU/Metal partitions**
for all 13,863 cells at each of three predeclared seeds. All twelve commands pass:
CPU and Metal publication plus source reconstruction for seeds 7, 19 and 42.
Both adjusted Rand index (ARI) and optimally matched assignment agreement are 1
in every CPU/Metal comparison, exceeding the predeclared 0.99 thresholds.
No seed or cell was excluded.

| Seed | CPU / Metal clusters | CPU-to-Metal ARI | Matched changed cells |
| ---: | ---: | ---: | ---: |
| 7 | 12 / 12 | 1.0 | 0 |
| 19 | 11 / 11 | 1.0 | 0 |
| 42 | 12 / 12 | 1.0 | 0 |

The source graphs are the previously qualified FP64 CPU and explicit FP32 Metal
outputs, with identical neighbor membership and edge coordinates but small
weight differences. Their identities and the executable were verified before
running. All plans retain resolution 1, maximum sweeps 100, maximum levels 32
and level tolerance 1e-7. The protocol and agreement gates were recorded before
clustering. Every original cell identity is checked against the PCA metadata.
Independent NumPy modularity reconstruction differs from the native values by
at most 1.67e-14. All six fits terminate with no aggregate moves and no
disconnected communities.

## Limits that remain

Seed choice still matters: within either backend, ARI is 0.82044 between seeds
7/19, 0.69750 between 7/42 and 0.65261 between 19/42. The backend comparison does
not prove a unique or seed-stable partition. No authoritative cell-type labels
are available in this Hagai qualification, so agreement with CPU is not evidence
of biological accuracy. Embedding behavior, other cohorts, cell-type/condition
preservation and broader performance remain separate gates. CPU FP64 stays the
default; this result does not promote Metal to a universal replacement.

## Evidence and reproduction

`run.py` runs the actual CLI and reconstruction, using the recorded source graph
and executable paths. `check.py` checks complete identities and cluster sizes,
recomputes modularity, evaluates every CPU/Metal comparison and retains all
within-backend seed comparisons. No new reference labels or graph weights are
inferred. `verify.py` verifies the retained archive and independently checks the
partition equivalence from all saved label vectors without rerunning the GPU.

The archive contains scripts, pre-run protocol/plans, all six clustering results,
receipts, logs and numerical checks. Large parent graph/source snapshots remain
externally hash-bound under `/Users/n/numivivo-metal-clustering-20260912`.
All 12 successful commands and their elapsed times are retained; this fixed-order
single execution is not a comparative performance benchmark.
