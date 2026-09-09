# Sparse multilevel Louvain clustering

Enable clustering alongside the existing reduction and neighbors options:

```json
"clustering": {
  "resolution": 1,
  "seed": 7,
  "maximumSweeps": 100,
  "maximumLevels": 32,
  "levelTolerance": 0.0000001
}
```

The same count-analysis artifact contains the partition and reconstructs it on
verification. No graph or representation is inferred when the required plan
options are absent. Older plans omit clustering and retain their behavior.

This implementation optimizes undirected weighted modularity on the native
fuzzy connectivity graph using multilevel Louvain. Nodes start in individual
communities at each level. A seeded SplitMix64 permutation determines node visit
order; candidate community IDs are considered in sorted order. Positive moves
must exceed 1e-12 in modularity. Sweeps continue until no node moves, then the
graph is aggregated, retaining internal edge mass on its diagonal. Original-cell
labels propagate through every level; final cluster IDs follow the first member
in retained source-cell order. No dense adjacency matrix is constructed.

Modularity is recomputed on the original sparse graph at every level. Objective
decrease, exhausted sweep budget, or exhausted level budget causes rejection
before a clustering receipt is published. The explicit level tolerance stops
further aggregation when its objective improvement is sufficiently small.
Resolution controls the modularity penalty and changes the resulting partition.
Seed, options, cell identities, cluster sizes, per-level sweeps/moves/objective,
termination reason and disconnected-community count are retained in the report.

Louvain is a heuristic and can produce disconnected communities. The report
exposes this property rather than silently splitting or renaming communities.
This implementation is not Leiden refinement and does not claim globally optimal
modularity or authoritative cell types. Repeatability for an identical input
is separate from stability across seeds, resolutions or biological datasets.

`check_clustering.py` reconstructs the same graph in NetworkX, independently
computes the native partition's modularity and connectivity, and compares three
NetworkX Louvain seeds (7, 19, 41). Exact partition equality is not an oracle for
this heuristic. Every reference objective and adjusted Rand index is retained.
A predeclared 0.02 maximum objective deficit flags a material reference gap; it
is an engineering comparison threshold, not biological acceptance. Declared donor
and condition associations are reported descriptively when multiple complete
levels exist. They neither prove batch correction nor establish cell types.

```sh
python check_cli.py --clustering --binary /path/to/numivivo --imported /path/to/imported --out /new/native
python check_clustering.py --report /new/native/report.json --out /new/reference.json
```

The real-data scope is the complete previously qualified Kang B-cell cohort
(2,651 cells) and PBMC3k (2,700 cells), using the default 20-component PCA and
15-neighbor graph. Count provenance, graph checks and resident/pair-work bounds
remain those of the existing benchmark and reduction owners.
[Native UMAP coordinates](EMBEDDING.md) now have separate reference checks.
Stronger partition stability/biological evaluation, Leiden refinement, donor-aware
integration and out-of-core execution remain open.

The seeded visit order is tied to the retained graph row order. Permuting rows
can lead to another local optimum; cell identities must be used when comparing
partitions or attaching labels to an external AnnData object.

Production evidence is under `evidence/2026-09-09-clustering`: 29 Swift tests in
seven suites, eight workflow assertions per dataset, and 18 existing CLI
assertions across 24 commands. Kang's production partition has eight communities
and modularity 0.619669; PBMC3k has eleven and modularity 0.754976. Independent
objectives agree within 6e-15, and neither partition has disconnected communities.
Reference ARI ranges are 0.563–0.670 for Kang and 0.851–0.882 for PBMC3k; the
Kang result is not partition-equivalent to the reference. These are descriptive
results, not biological or production promotion.

`publication-scope.json` identifies unrelated upstream genomic additions preserved
after the successful release check. The tested full binary is based on
`97c6f03` plus the fingerprinted clustering sources. This qualification does not
claim a later full-product rebuild including those upstream additions. Raw build
logs retain existing compiler warnings and their original whitespace.
