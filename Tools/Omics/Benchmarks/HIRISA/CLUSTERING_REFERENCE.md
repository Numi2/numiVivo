# Complete clustering references and runtime-aware verification

The buffered native full-cohort run subsequently completed publication, replay
and all independent checks: [complete result](FULL_CLUSTERING_RESULTS.md). The
original reference partitions and the evidence below remain unchanged.

The independent weighted Louvain references cover all **1,612,594 HIRISA cells**
and **34,707,084 directed symmetric graph edges**. Their original igraph 1.0.0
partitions and receipts remain unchanged. A separate NumPy calculation now
checks every saved label and includes every edge in each partition's modularity.
Full native clustering publication/replay remains a separate running gate at
this checkpoint; reference completion does not establish native completion.

| Frozen reference seed | Clusters | Saved modularity | Independent discrepancy |
| --- | ---: | ---: | ---: |
| 7 | 28 | 0.8875515259624742 | 6.695e-14 |
| 19 | 29 | 0.8864900168202439 | 8.904e-14 |
| 41 | 29 | 0.8830024043026399 | 7.749e-14 |

The original protocol uses resolution 1 and native seed 7, with maximum
100 sweeps, 32 levels, tolerance 1e-7 and 20 billion edge visits. Native modularity
may be at most 0.02 below the best of the three fixed reference partitions.
References use separate heuristic visitation/stopping rules, so identical
partitions are not required. Adjusted Rand comparisons retain that variability.
Neither a high modularity nor association with author labels validates cell types.

## Reusing numerical references across native runtimes

The optimized runtime must publish fresh PCA, graph and cluster receipts under
its actual identity. `reference_clustering.py check --reference-graph ...` can
apply the old numerical reference only when all seven graph payloads and all
six fitted PCA payloads match exactly in byte count and SHA256. Count-source
fingerprints must agree, each receipt must bind its actual payloads, and each
graph/PCA implementation pair must agree. Only fitted graph inputs are admitted.
The original graph-check receipt and igraph reference bindings remain unchanged.

The check validates all original cell identities, every label, cluster sizes
and canonical ordering, modularity over every edge, within-community
connectivity, level diagnostics and comparisons with all three reference seeds.
Before/after hashes reject changing input artifacts. This establishes the
applicability of the numerical references; native source reconstruction/replay
and biological acceptance are separate requirements.

`collect_clustering.py` hashes all nineteen required remote files before and
after collection. It transfers only the four cluster outputs and two new parent
receipts. The other thirteen files are hard-linked from an exactly matching
local qualified graph after checking each against the remote digest. It checks
the resulting bundle and runs the independent verifier. It neither copies raw
count sources nor fabricates receipt identities. A result marked collected or
independently passed does not assert that native replay completed.

## Qualification on complete real data

Complete Baron direct-reference and explicit cross-runtime checks pass for all
8,569 cells and produce identical partition metrics: 19 clusters, modularity
0.8728444955018703 and no disconnected communities. The maximum reference
modularity deficit is 0.0009669760858693754. Seven expected rejections cover
implicit cross-runtime reuse, rehashed graph/PCA changes, changed count-source
identity, inconsistent PCA implementation, unfitted input and changed reference
partitions. Original artifacts remain byte-identical after these checks.

Actual Mac mini → local collection transfers **556,051 logical payload bytes**
and reuses **17,069,386 bytes**; two remote snapshots and local hashes agree.
Attempting to collect Hagai with Baron parent data fails before creating a local
bundle or transferring its result. These are complete smaller-cohort collector
and verifier checks, not a completed optimized HIRISA partition.

The first test harness expected a content-mismatch diagnostic, but the rehashed
graph change correctly failed earlier at exact-size comparison. The failed
harness log and source remain retained. The corrected diagnostic expectation
passes without changing verifier behavior or numerical results.

The full-HIRISA NumPy reference check independently reconstructs all three
saved objectives within 8.904e-14. It reads edges in blocks and retains only
cell-scale labels/degrees, rather than allocating cells × genes. It does not
select replacement partitions or tune the original acceptance threshold.

## Reproduce and retain evidence

Use Python 3.11 or newer and the pinned
[clustering requirements](requirements-clustering.txt). Restore source-bound
graph and Baron bundles from the earlier [graph](evidence/2026-09-10-graph) and
[storage-access](evidence/2026-09-10-storage-access) archives. The current archive
records the exact external files and their upstream restoration locations for
both Baron runtime identities.

```sh
python check_reference_reuse.py --root /absolute/path/hirisa --out /new/checker-checks
python verify_clustering_reference.py \
  --graph /absolute/path/hirisa/graph/full \
  --protocol /absolute/path/hirisa/clustering/protocol.json \
  --graph-check /absolute/path/hirisa/graph/independent-check/checks.json \
  --reference /absolute/path/hirisa/clustering/independent-reference \
  --out /new/full-reference-check
```

For a completed new native publication with exactly matching parent payloads:

```sh
python collect_clustering.py --host macmini \
  --remote-bundle /absolute/remote/path/optimized-clustering/clustering \
  --reference-graph /absolute/local/path/graph/full \
  --protocol /absolute/local/path/clustering/protocol.json \
  --graph-check /absolute/local/path/graph/independent-check/checks.json \
  --reference /absolute/local/path/clustering/independent-reference \
  --out /new/collected
```

The [reference archive](evidence/2026-09-10-clustering-reference) contains all
three original full-HIRISA partitions, the independent full-edge check, Baron
checks, corruption logs, actual remote collection records and exact scripts.
Its 62 members use 2,232,485 compressed bytes; manifest SHA256 is
`d019ac29fade56486b5e0ff760e67877997886fc1aff2f4017d8a394de92443c`.
Verify with `python verify_archive.py evidence/2026-09-10-clustering-reference`.
Rebuild into a fresh directory with `archive_clustering_reference.py --root
/absolute/path/hirisa --output /new/archive`. The archive includes the collector
coordinator source; its presence does not assert a completed full native run.
