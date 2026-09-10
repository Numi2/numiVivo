# Source-bound cells for streamed aggregation

`singlecell-h5ad-pseudobulk` optionally accepts `cellSelection` in its plan:

```json
{
  "cellSelection": {
    "source": {"bytes": ["the 32 integer bytes of the source SHA-256"]},
    "observationIndices": [7, 2, 9],
    "provenance": "The explicit cohort rule and its evidence"
  }
}
```

The source fingerprint must match the copied H5AD snapshot. Indices are unique,
nonnegative, in bounds, and in the desired report order. A missing selection
keeps the existing complete-source behavior. An explicit selection must contain
at least one observation. The 2 MiB complete-plan limit still applies.

Swift/HDF5 scans and validates the complete count source, then routes selected
rows into its existing raw-count aggregation and per-cell quality accumulators.
No filtered H5AD or dense cells-by-features matrix is constructed. Changing legal
counts in excluded cells leaves the selected report unchanged. Invalid count
values in excluded cells are still rejected by source validation.

The bundle retains the entire original H5AD. Report `sourceObservationIndices`
maps report-local rows back to original rows, and `sourceCellCount` records the
original axis length. `metadata.cells`, `quality`, group `sourceCellIndices`,
and `canonicalNonzeros` describe the selected rows. All feature identities and
the original sample dictionary remain present, including unused samples. These
two new report fields are absent for a plan without selection, preserving the
previous default report format. CLI replay output distinguishes `sourceCells`
from `selectedCells` when selection is active.

The selection supports raw aggregation and its existing contrasts. Combining it
with per-cell program scoring or reduction currently fails explicitly because
those separate source passes still address the original row order. Use the
[axis projection](Projection/README.md) route when a materialized, aligned subset
is required and its documented input/output bounds admit the source.

`check_cell_selection.py` exercises CSR, CSC and dense input, nonmonotonic row
selection, zero-count retained cells, exact raw sums and mitochondrial metrics,
unchanged default behavior, excluded-count mutation, full-source validation,
source fingerprint binding, invalid/duplicate/empty indices, unknown fields,
unsupported reduction combinations, and receipt tampering. Its 26 checks pass
with the full release CLI. The 25 existing `check_streaming.py` checks also pass.
These fixtures establish software behavior; the separate
[Adamson cohort](../PerturbationPrediction/Adamson/COHORT.md) supplies experimental
count evidence.
