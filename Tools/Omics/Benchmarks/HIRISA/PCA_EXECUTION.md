# Complete HIRISA PCA execution specification

2026-09-10 at source commit 5cf4db930719628b123659d0e7426a13ab860850,
before HIRISA HVG selection or PCA fitting. Paired DE and donor-response results
have already been examined. This is a complete-source exploratory reduction,
not train-only preprocessing for a held-out prediction experiment.

Retain all 1,612,594 original cell identities, all 18,082 source genes and every
integer count from the verified H5AD with SHA256
0873e698ebf8770a54e6dba09724ffbeda5e1a67dbf24d4e223552d4aca7969c.
Use the already qualified native source mapping. Do not filter author labels,
subsample cells, pool donors or choose a response signature from DE outcomes.

Freeze the existing Seurat-style HVG/PCA settings: normalize each cell by its
own complete library to 10,000; log1p transformation; 2,000 requested variable
features; 20 mean bins; 20 components; maximum Krylov basis 128; seed 7;
relative residual tolerance 1e-6. Preserve projection centers and original
gene order, including the native rule for HVG ties. Do not lower feature or
component counts or weaken residual gates to admit this source.

First run the unchanged released native PCA command with its current maximum
2,000,000,000-byte cache and 20,000,000,000-entry-visit allowances. Retain any
typed source, cache, work, row-count or publication rejection. These are storage
and work admission limits, not scientific model settings. Independently measure
all-feature moments and the complete selected nonzero count before setting
larger resource allowances. Any repair must live in the shared native owner,
retain bounded windows and checked arithmetic, and requalify smaller complete
datasets plus the full HIRISA source. No dense cells-by-genes array is allowed.

Compare normalization and HVG statistics against independent chunked source
calculations and the pinned Scanpy selection behavior. Use a complete sparse
or file-backed linear operator for the independent centered PCA reference;
retain all source cells and normalization denominators. Check selected genes,
centers, variances, subspace agreement, native residuals and orthogonality,
and original row identities. Freeze output hashes before comparison. Publish
and replay the native result and record source passes, cache bytes, work,
elapsed time and peak resident memory. Keep software admission, full-source
numerical qualification and CPU/scverse performance comparison separate.

Graph, clustering and integration need separate frozen numerical settings and
biological-preservation comparisons after PCA qualifies. Author cell labels
remain predictions. This PCA experiment alone does not qualify annotation,
integration, unseen perturbations, Metal acceleration or cross-scale biology.

## Downstream admission audit during execution

The PCA storage repair does not yet admit complete HIRISA HNSW execution.
`Sources/NumiVivoCore/OmicsHNSW.cpp` independently rejects more than 1,000,000
rows and more than 128,000,000 directed neighbor entries. The Swift score/storage
readers now accept the larger axis, but that does not override the C++ gate.
Before the full graph stage, reconcile both owners, preserve overflow checks,
freeze graph/search settings, and measure recall against exact all-candidate
searches for a fixed query panel. No graph or integration completion follows
from the current PCA admission checks.
