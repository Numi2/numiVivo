# Native accessibility TF-IDF

`multiassay-10x-tfidf` normalizes an explicitly mapped accessibility assay while
retaining its count units, genomic intervals, feature namespace, observation
mapping and complete original source. RNA and unknown assay selections fail.

```sh
numivivo-omics multiassay-10x-tfidf original.h5 --plan mapping.json --assay atac --output new-bundle
```

Use the existing [paired Multiome mapping](../pbmc3k-multiome-plan.json). Build
the scoped actual CLI with `Tools/Omics/H5AD/build.sh <build-directory> --with-cli`.
This command uses the existing native 10x reader and a two-pass sparse kernel.
The complete multi-assay CSR input remains resident. Normalization retains only
cell/feature totals, ordering state and bounded digest/output buffers; it never
allocates cells by features. This is not out-of-core ingestion or Metal execution.

The fixed method is Signac method 1, scale 10,000:
`log1p((count / cellTotal) * (numberOfCells / featureTotal) * 10000)`.
The denominator uses total counts per feature, not the number of detected cells.
Counts are not binarized. Measured zero cells/features stay on their axes and
produce no sparse entries; missing modalities are not padded with measured zeros.
Both scans must supply identical ordered canonical entries. SHA256 comparison
rejects changes even when row and column totals remain unchanged.

Output contains `original.h5`, `mapping.json`, `axes.json`, `statistics.json`,
`values.bin` and `receipt.json`. Each sparse output record is little-endian
UInt32 assay-local row, UInt32 feature and Float64 value. Axes bind local rows
to source observations. The receipt binds all artifacts and implementation.
Publication stays private until the complete calculation succeeds; existing
destinations are preserved. The `multiassay-tfidf-verify <tfidf-bundle>` command now reconstructs the
complete normalization from the retained HDF5 source and checks every artifact.
The LSI command requires this verification and scans a private immutable values
snapshot. A deliberately changed value with a recomputed receipt hash is rejected
by source reconstruction; producer and verifier fingerprints remain distinct.

## Complete measured benchmark

The [paired 10x benchmark](../MULTIOME.md) supplies all 2,711 nuclei, 98,319 peaks,
19,292,713 nonzero measurements and 48,245,242 cut sites. There is one zero-total
peak and no zero-total cell. No cell or peak selection is applied.

- Native normalization matches the independent SciPy formula at every nonzero.
- The exact parsed `RunTFIDF.default` from Signac 1.16.0 executes in R/Matrix on
  the full matrix; maximum difference from the same oracle is 8.89e-16. This
  executes the upstream function, not the complete Signac/Seurat pipeline.
- The actual scoped product CLI produces byte-identical normalized records to
  the qualified native kernel. All totals, barcodes, feature IDs and receipt
  hashes match. RNA, missing-assay and existing-output cases are rejected.
- Fourteen kernel failure/cancellation checks pass, including duplicate and
  descending coordinates, overflows, invalid scales and altered replay with
  unchanged margins. Zero axes are preserved.

The source contains a single healthy donor. These results establish numerical
normalization and provenance, not LSI, peak calling, regulatory links, cell types,
joint RNA/ATAC integration, perturbation response or biological mechanism.
Peak memory and comparative end-to-end performance remain unqualified.

## Reproduction and retained evidence

[Evidence](evidence.tar.gz) and [manifest](manifest.json) retain scripts, kernels,
the qualified binaries, axes, statistics, receipts, checks and source hashes.
Every archive member was independently read and hash-verified. The 309 MB
normalized record file is hash-bound and retained in the local product bundle;
it is not embedded in the archive. The original 38,844,318-byte HDF5 is already
published with the paired benchmark and remains an explicit replay dependency.

Unpack into a new directory and adjust the recorded source, mapping and binary
paths in `check.py`, `check_signac.py` and `check_product.py` for that directory.
Run with Python/NumPy/SciPy/h5py, native HDF5 and R/Matrix available. `check.py`
streams both passes to the native harness and checks every output coordinate
and value without storing a dense matrix. `check_signac.py` executes the pinned
upstream function; `check_product.py` invokes the real command and verifies its
complete bundle. Preserve existing output directories and create a fresh one
for each run.

The initial R harness wrote to a literal file named `stdout`; its verified output
was retained and the pipe route corrected to `/dev/stdout`. The first product
build failed with ENOSPC; its log is retained alongside the successful retry.
Signac source is MIT licensed, copyright Tim Stuart; its notice and MIT template
are included. The 10x source is CC BY 4.0 and is not relicensed by this project.

Primary method: [Signac 1.16.0 source](https://raw.githubusercontent.com/stuart-lab/signac/1.16.0/R/preprocessing.R).
