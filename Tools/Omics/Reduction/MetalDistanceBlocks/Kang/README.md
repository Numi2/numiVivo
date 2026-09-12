# Full Kang Metal backend preservation and CLI timing

The explicit Metal FP32 neighbor backend passes the predeclared preservation
comparison on all **24,673 original Kang cells**, eight donor identities and
eight source-annotated cell types. Fresh PCA scores and metadata match the
original frozen cohort exactly. All fifteen native publication/reconstruction
commands pass. This uses fitted PCA, not the separately evaluated MNN integration.

CPU and Metal neighbor identities and order are identical. All 706,016 fuzzy
edge coordinates match; maximum common-edge weight difference is 1.304e-5.
At seeds 7, 19 and 42, clustering partitions agree exactly (ARI and matched
assignment agreement 1, zero changed cells), with 21, 20 and 21 clusters.
All source-type and condition-concordance regression gates pass without loss.
The pre-run thresholds were ARI/agreement ≥0.99 and loss ≤0.01 for source-type
ARI, each type's neighbor recall and stimulation balanced accuracy. Every type
has both source condition classes and remains in evaluation.

## Source-annotation concordance

The following values are identical for CPU and Metal. Votes use the 19 non-self
neighbors and fixed annotation-index tie breaking. These are local graph
concordance measures, **not held-donor prediction accuracy or causal effects**.

| Source type | Cells | Neighbor type recall | Condition balanced accuracy |
| --- | ---: | ---: | ---: |
| B cells | 2,651 | 0.987929 | 0.993975 |
| CD14+ Monocytes | 5,697 | 0.994032 | 0.996858 |
| CD4 T cells | 11,238 | 0.989678 | 0.978956 |
| CD8 T cells | 1,621 | 0.881555 | 0.959276 |
| Dendritic cells | 529 | 0.922495 | 0.998062 |
| FCGR3A+ Monocytes | 1,089 | 0.970615 | 0.991512 |
| Megakaryocytes | 132 | 0.992424 | 0.651829 |
| NK cells | 1,716 | 0.941725 | 0.968543 |

Overall condition balanced accuracy is 0.981816 in both paths. Cluster-to-source
type ARI is only 0.331232, 0.352184 and 0.313756 across the seeds; preserving it
does not validate automatic annotation. Megakaryocyte condition concordance
remains weak. The eight donor IDs are retained, but this evaluation does not
train a donor-held-out classifier or establish unseen-context generalization.
Source annotations never enter neighbor distances or clustering. CPU FP64
remains the default; no broader biological promotion follows from this result.

## Repeated complete CLI timing

After preservation checks, a fixed CPU/Metal/Metal/CPU/CPU/Metal sequence ran six
fresh native graph publications. All eight graph/plan/execution/receipt files
match the corresponding qualified backend in every repetition.

| Backend | Process elapsed seconds | Median |
| --- | --- | ---: |
| CPU FP64 | 9.958, 10.014, 10.074 | 10.014 s |
| Metal FP32 | 6.968, 6.931, 7.018 | 6.968 s |

Metal median elapsed time is **30.4% lower (1.44x)** for this command and cohort.
Timing includes CLI startup, parent PCA source reconstruction, graph construction,
hashing and publication. Post-run comparison hashing is outside the interval.
CPU peak RSS ranges 169.0–203.3 MB; Metal 179.1–206.6 MB. Thus no memory reduction
is claimed. Caches were exercised on a shared physical M4 Pro desktop. All runs
are retained; fixed order and three repetitions do not establish a confidence
interval, cold-cache performance or universal speedup. No Scanpy performance
comparison, million-cell execution or acceleration of all single-cell stages is
claimed. These results extend the earlier Hagai research and native checks.

## Evidence

`run.py` freezes source identities, configuration and preservation gates before
actual CLI execution. `check.py` compares every cell, neighbor, edge coordinate,
partition and source-label metric. `benchmark.py` records its sequence and
identities before repeated execution and hashes every output. The archive
retains all plans, original result vectors, receipts, logs and comparisons.
Parent sources/graphs and binary arrays remain externally hash-bound under
`/Users/n/numivivo-metal-kang-20260912`. `verify.py` checks retained evidence and
partition equivalence; it does not rerun prediction, GPU work or scoring.
