# Sparse ATAC LSI development qualification

The native LSI kernel reuses the existing residual-qualified Krylov solver with
zero centering to decompose the complete TF-IDF matrix. It returns singular
values, left singular vectors, standardized cell embeddings and feature loadings.
The normalization denominator is unchanged. All peaks and all components remain
present; no depth-correlated component is silently removed.

The retained complete paired 10x dataset has 2,711 nuclei and 98,319 peaks.
The development protocol fixed 30 components, a 256-vector basis, seed 7 and
relative eigen-residual tolerance 1e-5 before execution. The actual native
library was linked into a scoped harness; no solver stand-ins were used.
All 30 components passed comparison against independent SciPy ARPACK SVD:

| Check | Maximum difference |
| --- | ---: |
| Relative singular value | 4.78e-15 |
| Left singular vector, after sign alignment | 1.79e-9 |
| Feature loading, after sign alignment | 9.29e-11 |
| Standardized embedding, after sign alignment | 9.32e-8 |
| Independent right singular-triplet residual | 2.08e-10 |

The components retain 14.12% of uncentered matrix energy. This is not centered
explained variance. The first component correlates with log1p total cut sites
at r=0.9629; it must not automatically be interpreted as biological variation.
The entire source contains one donor and has no perturbation validation.

Execution makes 602 operator scans plus the initial validation/energy scan.
Sparse source bytes are mapped in the harness; basis vectors and small output
axes remain resident. This is numerical qualification, not a comparative
performance result or million-cell memory qualification.

## Product execution and evidence

Build the actual scoped CLI using `Tools/Omics/H5AD/build.sh <build> --with-cli`.
Create a plan with all five fields:

```json
{"schemaVersion":1,"components":30,"maximumBasis":256,"relativeResidualTolerance":0.00001,"seed":7}
```

```sh
numivivo-omics multiassay-tfidf-verify <tfidf-bundle>
numivivo-omics multiassay-lsi <tfidf-bundle> --plan <plan.json> --output <new-bundle>
```

The complete product run passed all output values, matrix coordinates, source
axes and seven artifact fingerprints. Singular values, U and loadings match the
initial native result exactly; standardized embeddings agree within 1e-12
relative and absolute tolerance after the stable variance repair below.
The bundle contains the input receipt, axes, plan, model diagnostics, U,
standardized embeddings, feature loadings and its own receipt. Binary matrices
contain every coordinate as little-endian UInt32 row, UInt32 column and Float64
value. Model diagnostics include per-component depth correlations.

Input TF-IDF verification reconstructs normalization from the retained HDF5 and
uses a private snapshot of the value stream. Unknown plan keys, wrong schema,
invalid basis/components, self-rehashed TF-IDF tampering and existing output
destinations all passed rejection/preservation checks. Output is published only
after successful completion. A separate LSI output replay command is not yet
implemented.

An edge test exposed roundoff in the original naive variance calculation:
constant components could receive nonzero standardized values. Welford variance
now rejects these components. All 16 constant-matrix cases (2 through 17 cells)
passed against the final product library. The initial source/result and corrected
source remain distinct in the evidence.

`evidence.tar.gz` and `manifest.json` retain the frozen protocol, initial native
result, independent SciPy result, comparison and admission scripts, product
outputs, source/build fingerprints and logs. Original HDF5 and full TF-IDF values
are identified by their existing upstream evidence; they are not duplicated in
this archive. Scripts record the development host paths and require those inputs
or equivalent restored copies for replay.

This milestone does not establish joint RNA/ATAC integration, held-out biological
prediction, million-cell scaling or clinical validity.

The embedding convention follows the [Signac 1.16.0 source](https://raw.githubusercontent.com/stuart-lab/signac/1.16.0/R/dimension_reduction.R):
unit left singular vectors, optionally standardized per component. The native
result exposes both forms without clipping or component exclusion.
