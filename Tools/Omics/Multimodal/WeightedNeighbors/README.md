# Native paired weighted neighborhoods

The `VivoWeightedNeighbors.fit` API connects native cell-specific modality
weights to exact, directed weighted-neighbor selection and a UMAP-compatible
fuzzy-union graph. Both input representations retain their independent feature
spaces and exact paired cell identities.

All 2,711 measured nuclei use the unchanged 30-component RNA and ATAC inputs
from [PairedNeighbors](../PairedNeighbors/README.md). The weighting stage uses
the [qualified SNN-bandwidth calculation](../ModalityWeights/README.md).
Its existing exact neighbor results are reused during weighted selection.

For each ordered cell pair, each modality contributes its weight times a
distance-decaying affinity, after subtracting that cell's nearest-neighbor
distance. The distance is the square root of half the weighted dissimilarity.
The kernel power is one. Every other cell is considered; 14 non-self neighbors
are selected and self occupies slot zero, preserving the existing 15-slot
UMAP graph convention. Ties use ascending cell index.

This follows the affinity calculation in
[Seurat 5.3.0 MultiModalNN](https://github.com/satijalab/seurat/blob/v5.3.0/R/clustering.R).
It deliberately reports exact all-candidate search, rather than claiming parity
with Seurat's approximate candidate union. Its final fuzzy-union graph is also
distinct from Seurat's shared-neighbor graph. No full Seurat execution is claimed.

## Numerical qualification and repaired failure

Independent SciPy distance blocks and previously independently calculated
modality weights check every selected neighbor and distance. UMAP-learn 0.5.12
checks every fuzzy graph weight.

| Check | Corrected result |
| --- | ---: |
| Neighbor-index mismatches | 0 |
| Maximum distance error | 3.61e-16 |
| Maximum fuzzy-weight error | 2.68e-6 |
| Directed CSR entries (symmetric graph) | 55,398 |
| Connected components | 1 |
| Isolated cells | 0 |

The initial implementation failed despite exact neighbor identities: subtracting
affinity from one caused cancellation near zero distance. Maximum distance error
was 1.054e-8 and final graph-weight error was 0.7373. The corrected kernel computes
the equivalent weighted `-expm1(-distance/bandwidth)` expression. The independent
checker also derives nearest distances from its own distance blocks, avoiding a
mixture of distances from different floating-point implementations. The original
failed source, outputs and comparison remain retained; acceptance thresholds
remain 1e-10 for distances and 1e-5 for graph weights.

Modality swapping leaves the complete native graph unchanged. The weighting
refactor retains exact previous product results and passes the previous six
invalid-input checks. The corrected product build is separately checked against
the qualified corrected result before publication.

## Structural comparison and limits

| Baseline | Mean shared fraction of 14 neighbors | Cells with no shared neighbors |
| --- | ---: | ---: |
| RNA-only | 29.90% | 375 |
| ATAC-only | 66.82% | 2 |
| Equal-weight joint | 60.05% | 13 |

These are structural comparisons, not biological-preservation scores. The result
is closer to ATAC neighborhoods by this measure; that does not establish
biological superiority. The depth-correlated first ATAC component remains
included. No labels are treated as authoritative and the single donor cannot
qualify donor integration or unseen-context prediction.

The implementation retains bounded neighbor heaps and sparse graph structures,
not a dense cell-pair matrix. Exact-search limits still apply; this does not
qualify million-cell performance.

The archive and manifest retain the protocol, source/build fingerprints,
complete corrected and failed outputs, independent checker, native harnesses
and regression evidence. Input paths refer to the prior published paired
studies and must be restored or adjusted for replay. Next work is a reproducible
product bundle/CLI boundary and biological-preservation evaluation against
separately specified experimental evidence. General perturbation prediction,
regulatory prediction and cross-scale biological qualification remain open.
