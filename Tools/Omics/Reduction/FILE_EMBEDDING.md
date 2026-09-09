# File-backed UMAP-compatible embedding

The binary PCA graph store feeds the same numerical optimizer as the resident
count-analysis route. Edges are scanned through 16 MiB windows. Each retained
edge's head, tail, period and next positive/negative sample times occupy a
32-byte record in private scratch storage. One shared 16 MiB mapping window
preserves schedule updates between sequential epochs without resident edge arrays.

```sh
numivivo singlecell-graph-embed /path/to/binary-graph \
  --plan embedding-plan.json --output /new/embedding-bundle
numivivo singlecell-graph-embed-verify /new/embedding-bundle
```

```json
{
  "schemaVersion": 1,
  "embedding": {
    "dimensions": 2,
    "epochs": 500,
    "minimumDistance": 0.1,
    "spread": 1,
    "learningRate": 1,
    "negativeSampleRate": 5,
    "repulsionStrength": 1,
    "seed": 7,
    "maximumUpdates": 5000000000
  }
}
```

Fitted and frozen-query inputs from exact or HNSW graph construction are
supported. The parent graph must use `storage: "binary"`. Publication retains and
reconstructs the complete graph/PCA lineage, binds result and execution hashes
to the executing binary, and atomically publishes the new bundle. A new binary
must refit inputs from original sources. Verification replays the entire lineage
and optimizer; changed results cannot be made valid by rehashing receipts.

Frozen-query inputs fit a new layout of the query cells; they do not place cells
into an existing reference UMAP. Real-cohort quality qualification covers fitted
2D layouts. Query and 3D paths have lifecycle and trajectory regressions, without
corresponding real-data quality claims.

Only the requested two or three PCA columns are retained for initialization.
All score coordinates are validated in 2,048-row tiles. Cell identities,
initial/final coordinates, initialization columns and result JSON remain resident.
This removes edge-scale schedule storage from memory, not every cell-scale
allocation. No cells-by-genes matrix is built.

Curve fitting, PCA scaling/jitter, edge threshold/order, floating-point schedule
increments, clipped attraction/repulsion, learning-rate schedule and SplitMix64
sampling retain the existing formulas. Scratch is not a restart checkpoint and
is removed before publication and on failure. The [original documentation](EMBEDDING.md)
describes the numerical model and interpretation limits.

The default `maximumUpdates` remains 200 million. Larger runs explicitly request
a higher budget, now capped at 20 billion. The cap bounds both retained-edge
visits and conservative/actual attractive-plus-negative updates. Insufficient
work rejects the request; epochs are never silently reduced. Input PCA/graph
work budgets remain separate. Scratch requires 32 bytes per retained directed
edge, in addition to retained lineage and temporary replay disk space.

## Qualification design

`check_file_embedding.py` covers fitted/query and exact/HNSW inputs, two and three
dimensions, byte-identical replay, frozen numerical-oracle equality, rehashed
coordinates/execution/parent rejection, low work budgets, JSON parents and overwrite.
Native regressions force multiple windows and check complete trajectories and
cleanup on callback failure.

`LegacyEmbeddingReference.swift` freezes the result shape and numerical algorithm
from `36ccd85`; it is built separately and never linked into the product. It
imports current public input models, including expanded budget admission.
Equality proves preservation of the numerical algorithm, not that the old CLI
accepted these larger-budget plans.

The predeclared protocol uses 500 epochs and two dimensions for full Baron,
Hagai and Norman cohorts. Every input graph coordinate/FP64 value, PCA score byte
and cell identity must equal the previous independently qualified HNSW graph.
Complete result JSON must equal the frozen numerical oracle. Independent checks
verify curve fitting and all sample counts, then run umap-learn's Float32/Tausworthe
optimizer with seeds 7, 19 and 41 on the same graph and initial coordinates.

Neighborhood quality is evaluated at 15 neighbors. Every Baron/Hagai cell is
checked; Norman uses 2,048 fixed, preselected queries against all 111,445 cells.
The same queries score initialization, native output and all references. Distance
calculations use 64-row blocks, avoiding a full dense cells-by-cells matrix.
Gates require native trustworthiness at least 0.90, within 0.02 of the lowest
reference trustworthiness, and recall within 0.05 of the lowest reference recall.
These engineering gates do not establish biological preservation, optimizer
convergence, meaningful density or authoritative cell types.

Million-cell execution, file-backed integration, transformation of unseen cells,
spectral initialization, biological/cross-donor evaluation and Metal/scverse
end-to-end comparisons remain open. Fixed-epoch coordinates are descriptive;
clustering continues to use the PCA graph.

## Full-cohort measurements

The [evidence archive](evidence/2026-09-09-file-embedding/README.md) retains
complete native trajectories, exact frozen-oracle results, every input graph/PCA
check, all three independent optimizer coordinate sets, fixed query rows, work
counts, source/binary hashes and raw command logs. Each cohort used 500 epochs.

| Cohort | Cells | Quality queries | Native trustworthiness | Reference range | Native 15-neighbor recall | Reference range |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Baron | 8569 | 8569 | 0.98343 | 0.98276–0.98307 | 0.28920 | 0.28492–0.28654 |
| Hagai | 13863 | 13863 | 0.96461 | 0.96465–0.96528 | 0.13929 | 0.14034–0.14081 |
| Norman | 111445 | 2048 | 0.91497 | 0.91537–0.91633 | 0.03988 | 0.03991–0.04121 |

All predeclared quality gates pass. Norman values estimate quality from the
fixed query sample; they are not exhaustive all-cell rank measurements. Every
optimizer still used the full cohort and graph. Low 15-neighbor recall remains
visible: the two-dimensional layout loses most original local relationships.
No cell labels or donor-integration conclusions follow from these coordinates.

| Cohort | Publish seconds | Publish peak resident bytes | Verify seconds | Verify peak resident bytes |
| --- | ---: | ---: | ---: | ---: |
| Baron | 12.88 | 210239488 | 12.97 | 212156416 |
| Hagai | 25.38 | 231096320 | 25.53 | 228786176 |
| Norman | 254.56 | 486113280 | 256.77 | 482230272 |

[Complete native/reference layouts](evidence/2026-09-09-file-embedding/figure/coordinates.png) show all cells without inferred labels; source hashes are retained alongside the figure.

These are single-run CPU observations on M4 Pro/24 GiB/macOS 26.6. Each publish
and verify includes complete parent PCA/graph reconstruction and serialization.
The frozen numerical oracle starts from an existing graph and scores; its timing
is not a comparable workflow measurement. Independent references ran on the
laptop. Completed artifacts were copied while later native commands ran, without
a competing native build or benchmark. This is not a same-host scverse or GPU
speed comparison.

Both builds passed, along with 83 Swift tests in 25 SingleCell/MultiAssay/CountStore
suites, 12 full-cohort production commands, three numerical-oracle runs and 88 CLI
fixture commands including 31 expected rejections. The archive distinguishes the
new explicit work-budget admission from the preserved numerical algorithm.
