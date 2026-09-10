# Complete HIRISA graph qualification

The full native binary graph qualifies every one of the 1,612,594 original cells
under the frozen numerical and sampled-recall protocol. Publication, native
reconstruction, every-record independent numerical checks and both recall gates
pass. The complete PCA prerequisite is published as `e4b630b6`.

## Frozen source, algorithm and execution

The [execution protocol](GRAPH_EXECUTION.md) fixes all cells, 20 PCs, k=15,
connections=16, construction width=200, search width=128 and seed=7. The graph
uses the exact source PCA score bytes already independently qualified in
[PCA results](PCA_RESULTS.md). The new executing binary is
`ee332d69aa3b3514444f00c4474f03332c296e236ed61d26cf2e2fb933e3d699`;
its manifest binds all 532 authored production source files, compiler, HDF5 and
Mac mini hardware. No executable or receipt identity was relabeled.

A complete new-runtime PCA refit took 577.5855 seconds. Plan, model, metadata,
quality, scores and loadings all equal the previously qualified bytes exactly.
Native graph publication reconstructs that source again before graph construction.

| Observed native publication quantity | Value |
| --- | ---: |
| Cells / PCs / k including self | 1,612,594 / 20 / 15 |
| Neighbor records, including self | 24,188,910 |
| Directed symmetric connectivity records | 34,707,084 |
| Connected components / isolated cells | 1 / 0 |
| Construction metric evaluations | 4,208,690,951 |
| Query metric evaluations | 2,755,059,586 |
| Total metric evaluations | 6,963,750,537 |
| Explicit metric allowance | 16,000,000,000 |
| Actual decoded score cache bytes | 258,015,040 |
| Explicit score cache allowance | 258,048,000 |
| Encoded score read bytes / calls | 516,030,080 / 6,300 |
| Estimated serialized HNSW index bytes | 245,941,096 |
| Whole native publication seconds | 1,097.1219 |
| Maximum publication resident bytes | 4,360,978,432 |
| Whole native reconstruction seconds | 1,055.6142 |
| Maximum reconstruction resident bytes | 4,351,361,024 |

The score cache fits the complete compact PCA representation in this explicit
plan. The index and metadata remain resident. The index-size estimate excludes
allocator, mutex and container overhead; whole-command RSS includes source
reconstruction. Brief smaller-cohort clustering regressions also ran on this
host during source reconstruction. These timings are observations, not a
controlled CPU/scverse or GPU comparison.

## Completed prerequisite evidence

- Twelve native tests in four suites, including the real C++ million-row
  admission/work-budget test and cache eviction regression, pass.
- Complete release product builds successfully.
- Fitted/query HNSW lifecycle and binary graph-store lifecycle pass, including
  reconstruction, exact repeats and rehashed tampering controls.
- Every binary graph coordinate and FP64 value, source score byte and cell
  identity in complete Baron and Hagai agree with the earlier independently
  qualified graphs. Their native graph reconstructions pass.
- The independent full binary graph checker passes end to end on Baron,
  including umap-learn fuzzy topology, weights and every returned distance.
- The frozen 2,048-query panel was searched exactly against all 1,612,594
  candidates before the HIRISA graph was inspected: 3,302,592,512 FP64 SciPy
  distance comparisons. Preparation took 25.5972 seconds on the laptop; it is
  not equivalent to the native all-cell graph workload.

## Independent complete graph and frozen recall checks

`reference_graph.py` passes every returned distance, binary coordinate, source
score/metadata/QC byte identity, bandwidth, kernel residual and fuzzy edge check.
The independently reconstructed graph also has one component and no isolated
cells. Fuzzy topology matches umap-learn exactly.

| Independent measurement | Result |
| --- | ---: |
| Fixed exact-query cells | 2,048 |
| Candidates searched per query | 1,612,594 |
| Mean strict recall | 0.9997209821428572 |
| Fifth-percentile strict recall | 1.0 |
| Minimum strict recall | 12/14 |
| Mean tie-aware recall, reported separately | 0.9997907366071429 |
| Maximum returned-distance error | 1.7763568394002505e-15 |
| Maximum fuzzy-weight error | 5.189564978502759e-5 |
| Whole independent-check seconds | 10.8088 |
| Maximum independent-check resident bytes | 5,302,026,240 |

Both unchanged strict gates pass: mean ≥0.95 and fifth percentile ≥0.85. This
is sampled exact recall against the complete candidate universe, not exact
recall for every query cell. All 34,707,084 fuzzy edge records meet the frozen
mixed tolerance `absolute error ≤ 1e-5 + 1e-4 × abs(reference weight)`; the
maximum absolute error is not compared with the absolute term alone. Every
returned distance meets rtol/atol=1e-12. The independent sparse reference graph
is resident and runs on the laptop; this is not a paired performance comparison.

## Complete reconstruction and archive

Native full graph replay passes with the frozen executable and original OS
identity, including complete source reconstruction. The graph receipt SHA256 is
`6db0d3deb8ae66bc20108edb57676b8f24a0bddae8e938214d189d9b483e4073`.

The [complete archive](evidence/2026-09-10-graph/manifest.json) contains 691
members covering 649 logical source files. Its manifest SHA256 is
`79020ba5c5db5adef4467037d665308c8d99a0bfc1f64a3be7b5d9792383da10`.
Run `python verify_archive.py evidence/2026-09-10-graph` to independently verify
every encoded member, decoded chunk and reconstructed full-file hash. The archive
reuses verified compressed PCA chunks; only 714,334,830 new encoded bytes are
needed beyond the earlier PCA archive. Original large source H5ADs and executables
remain external with their immutable identities recorded.

Local transport reuses the six byte-identical PCA data files under the freshly
executed native input receipt. The original read-only transfer was stopped to
avoid sending duplicate gigabytes; the native computation was not stopped.
Transport receipts and all interrupted-transfer records are retained.

Clustering has started after successful graph reconstruction and independent
qualification; its full native and independent acceptance remains separate.
Integration still has separate algorithm and artifact scale gaps.
This graph does not establish biological signal preservation, authoritative cell
annotation, a fully out-of-core workflow, Metal performance or the broader goal.
