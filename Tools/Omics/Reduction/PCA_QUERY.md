# Frozen PCA projection of new cells

`singlecell-h5ad-pca-query` applies a [standalone PCA bundle](PCA_BUNDLE.md) to
new cells without selecting query HVGs, estimating query centers or refitting
loadings. The same `VivoFrozenPCAProjection` calculation now serves the existing
resident reference-label mapper. This command emits coordinates and QC only.

```sh
numivivo singlecell-h5ad-pca training.h5ad --plan fit.json --output training-pca
numivivo singlecell-h5ad-pca-query query.h5ad --plan query.json \
  --reference training-pca --output query-pca
numivivo singlecell-h5ad-pca-query-verify query-pca
```

The fit plan must declare a nonempty `featureNamespace`, such as
`HUMAN_GENE_SYMBOL`. Older standalone fits without a namespace remain valid for
PCA reconstruction but cannot serve as query references. Query plans contain
`schemaVersion: 1`, the explicit H5AD `mapping`, the same `featureNamespace`, and
optional `maximumProjectionUpdates` (default 2 billion, maximum 20 billion).
Query group/label mapping is rejected. Namespaces are user-declared provenance;
no automatic gene-symbol synonym resolution, identifier conversion or liftover
is performed.

The complete query feature-ID set must equal the training set, even for genes
outside the selected HVGs: every gene contributes to the normalization total.
Feature order may differ. Organism and count units must match. Overlap by exact
`sampleID`/`barcode` identity is rejected, and overlapping donor IDs are reported.
These checks do not discover renamed duplicate cells or certify an independently
authored experimental split. The caller must preserve source identities and
choose the training population before fitting.

Each query library is normalized using the training target, followed by log1p.
Selected values are multiplied by frozen loadings and shifted by training means:
`score = log1p(target * counts / all_gene_total)[selected] @ loadings - centers @ loadings`.
A zero-count cell retains the mathematical centered-zero score and is explicitly
flagged by its QC total and the report's `emptyLibraries`; it is not assigned a
biological label. No query mean is subtracted.

## Storage and reconstruction

The result retains `original.h5ad`, canonical `plan.json`, `metadata.json`,
`quality.json`, `report.json`, `scores.bin`, `receipt.json` and a complete immutable
`reference/` PCA bundle. Binary scores use the same complete row-major 16-byte
record format as standalone training PCA. The report records dimensions, selected
entries, scalar projection updates, scratch bytes, source passes, empty libraries
and donor overlap. The receipt binds every artifact and the reference receipt.

Publication first snapshots and reconstructs the training bundle. Query execution
then scans sparse HDF5 twice: the first pass obtains exact all-gene QC and checks
the projection work budget; the second accumulates selected-feature contributions.
Private row-major FP64 scratch maps at most 16 MiB in page-aligned, row-aligned
windows. Scores are emitted with a 1 MiB write buffer. Scratch is removed before
atomic publication. Query metadata, QC, feature lookup and selected loadings remain
resident; training reconstruction still has its own resident fitting costs.

Verification reconstructs the retained training fit and the entire query result.
Rehashing edited loadings, score coordinates/values or report work counts does not
bypass reconstruction. Verification includes training-fit cost and is not a
measurement of the query transform alone. Retaining a full training snapshot per
query uses disk space; no shared immutable artifact-store optimization is claimed.

## Evidence tools

- `check_pca_query.py`: CSR/CSC and reordered-feature fixtures, repeatability,
  zero libraries, overlap/namespace/units/organism/label/budget rejection and
  rehashed training/query tampering.
- `run_pca_query_baron.py`: all four complete Baron donor-held-out folds, native
  query reconstruction and existing reference-mapper regression.
- `pca_query_reference.py`: every binary coordinate, source identity, exact QC,
  independent Scanpy normalization and SciPy frozen projection, plus optional
  exact comparison with legacy native query scores.

The Python programs are external qualification tools. They materialize arrays;
the native query path does not. Passing projection checks establishes correct
application of a training fit, not donor integration, perturbation prediction,
biological generalization, million-cell throughput or Metal acceleration.

## Complete Baron qualification

The [archived qualification](evidence/2026-09-09-pca-query/README.md) covers every
Baron human cell as a held-out query, using all 20,125 genes and training-derived
2,000-HVG/20-PC models:

| Held-out donor | Query cells | Publish seconds | Publish maximum resident bytes | Verify seconds |
| --- | ---: | ---: | ---: | ---: |
| human1 | 1,937 | 5.25 | 170082304 | 5.25 |
| human2 | 1,724 | 5.33 | 156712960 | 5.32 |
| human3 | 3,605 | 4.84 | 167182336 | 4.85 |
| human4 | 1,303 | 5.41 | 165019648 | 5.47 |

All 8,569 query cells are disjoint from their training fold. Independent
Scanpy/SciPy score error is at most 2.85e-14 absolute. Training numerical results
and every legacy reference-mapping report remain exact. All 71 Swift tests in
20 suites, 20 real-data native commands, 22 fixture commands (13 expected
rejections) and the existing standalone PCA CLI regression pass. Publication and
verification times include training reconstruction; the legacy label mapper
performs additional classification work and is not an equivalent timing baseline.
