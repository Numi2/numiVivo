# Cell-specific modality weights

The native `VivoModalityWeights.fit` API now computes paired RNA/ATAC
embedding-predictability weights. All 2,711 measured nuclei pass independent
NumPy/SciPy comparison for bandwidths, reconstruction distances, scores and weights.
The product library reproduces the complete standalone result exactly.
Swapping modality inputs swaps weights exactly. Six invalid-input checks cover
row counts, non-unit rows, nonfinite values, duplicate identities, invalid k
and degenerate constant representations.

The input is two explicitly paired, finite, row-L2-normalized representations.
The current benchmark uses the unchanged 30-component inputs from
[PairedNeighbors](../PairedNeighbors/README.md), including the first ATAC component.
No cell selection, biological label supervision or feature-space merging occurs.

The algorithm follows the two-modality weighting calculations in
[Seurat 5.3.0 clustering.R](https://github.com/satijalab/seurat/blob/v5.3.0/R/clustering.R)
and its SNN-width helper. The frozen configuration uses 15 neighbors including
self, unsmoothed scores, shared-neighbor far-distance bandwidth, scale 1, cross
constant 1e-4 and clipping to [0,200]. Exact native neighbors replace approximate
search; full Seurat execution and approximate-search parity are not claimed.
The independent reference constructs sparse shared-neighbor intersections with
SciPy and computes all intermediates separately. Source links and hashes are
retained; upstream source files are not redistributed in this archive.

| Intermediate | Maximum absolute error |
| --- | ---: |
| SNN bandwidth | 4.45e-16 |
| Within-modality distance | 2.23e-16 |
| Cross-modality distance | 4.45e-16 |
| Clipped score | 2.14e-13 |
| Modality weight | 4.45e-16 |

All checks pass the frozen 1e-10 absolute tolerance. RNA weights have mean
0.3740 and median 0.4270, ranging from approximately 5.57e-87 to 0.9619.
These are embedding-predictability weights, not biological importance or
confidence in cell types. Some scores saturate at the declared clipping bound;
the full per-cell diagnostics remain available.

The implementation uses inverted neighbor membership for nonzero shared-neighbor
intersections, avoiding a dense cells-by-cells SNN. Existing exact search has a
bounded distance-pair budget. No million-cell or performance claim follows.

## Evidence and remaining integration

The archive contains the frozen input/source fingerprints, actual source,
standalone and product test binaries, independent reference, all numerical
results, test source and build logs. The actual scoped product library includes
the new owner via `Tools/Omics/H5AD/build.sh`; its complete source fingerprint
list and library fingerprint are recorded.

A local out-of-space write failure occurred before any native run. Verified
published transport duplicates and the archived initial LSI result were removed;
raw data and current paired inputs were retained. Development paths in scripts
must be restored or adjusted for replay.

Next is weighted neighbor selection and graph construction, compared against the
fixed RNA-only, ATAC-only and equal-weight baselines. A product bundle/CLI boundary,
biological-preservation qualification and general perturbation prediction remain
unfinished. These weights alone must not be presented as completed WNN integration.
