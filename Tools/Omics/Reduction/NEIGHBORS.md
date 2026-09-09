# Exact PCA neighbors and a UMAP-compatible fuzzy graph

Add `neighbors` alongside `reduction` in the existing analysis plan:

```json
"neighbors": {
  "neighbors": 15,
  "localConnectivity": 1,
  "maximumDistancePairs": 50000000
}
```

The optional `neighbors` result is carried through the same analysis artifact,
export and native reconstruction as PCA. It requires an explicit reduction plan;
there is no automatic change of representation or fallback to raw gene counts.
The default remains no graph for older plans.

The fixed-width neighbor arrays include self in slot zero. Remaining slots are
ordered by Euclidean distance, then retained cell index, giving deterministic
behavior for ties and duplicate score vectors. Each unordered cell pair is
measured once; bounded per-cell max heaps retain the exact k nearest entries.
There is no dense cell-by-cell distance matrix. Runtime is quadratic in cells;
this is an exact baseline, not an approximate or million-cell search algorithm.
The explicit pair budget is checked before distance computation. Its default
permits 50 million unordered pairs and its configurable upper bound is 500
million. Neither a performance claim nor an automatic budget increase is made.

Local fuzzy memberships use the UMAP reference convention: rho from positive
neighbor distances, a bandwidth search targeting log2(k), 64 search iterations,
1e-5 mass tolerance, bandwidth 1, and a 1e-3 mean-distance bandwidth floor.
Memberships are combined by fuzzy union `a + b - a*b` (set-operation mix ratio 1).
The result is symmetric connectivity CSR with sorted columns and no diagonal.
Weights, per-cell rho/sigma, post-floor kernel-mass residuals, connected component
count, isolated cell count, score dimensions and cell identities are retained.
Duplicate points or the bandwidth floor can prevent exact target mass; residuals
remain visible rather than being relabeled as convergence.

Native calculations use Double. `umap-learn` uses Float32 for parts of its kernel
calibration, so the reference check declares absolute 1e-5 / relative 1e-4 weight
tolerances. It requires exact neighbor membership/order and exact CSR structure,
then checks every edge weight, bandwidth, symmetry, diagonal and component count.
Independent exact distances use SciPy `cdist` in blocks of at most 64 cells;
reference code also avoids a complete dense distance matrix.

```sh
python check_cli.py --neighbors --binary /path/to/numivivo --imported /path/to/imported --out /new/native
python check_neighbors.py --report /new/native/report.json --out /new/reference.json
```

Real-data comparisons use all Kang B cells (2,651 cells; 20 PCA dimensions) and
PBMC3k (2,700 cells; 20 dimensions), with 15 neighbors including self. Original
count and PCA qualification remain in the benchmark and reduction evidence.
The graph tests are numerical evidence, not biological validation. Native UMAP
coordinate optimization, clustering, batch integration, approximate search and
out-of-core execution remain separate open requirements.

The reference algorithm is documented in the installed, version-pinned
`umap.umap_` functions `smooth_knn_dist`, `compute_membership_strengths` and
`fuzzy_simplicial_set`. Package versions are recorded in each comparison report.

Published qualification evidence is under `evidence/2026-09-09-neighbors`: 26
Swift tests in six suites, seven production assertions per dataset, and 18
existing CLI assertions across 24 commands. All neighbor identities/order and
sparse edge positions match the independent reference. The graphs have 56,420
(Kang) and 57,880 (PBMC3k) directed CSR entries, with maximum weight errors of
5.705e-6 and 3.885e-6. Both have one connected component. These are descriptive
numerical properties, not evidence of preserved donor or biological structure.
Raw compiler logs retain pre-existing warnings and their original whitespace.

The `integrated` evidence subdirectory records the final rebuilt executable
after preserving upstream genomic-document validation commit `70f5e1e`. The
parent evidence directory retains the earlier successful run. Both identify
their exact base commit, unchanged neighbor source hashes and executable hash.
