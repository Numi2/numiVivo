# Adamson UPR: complete native ingestion, prediction qualification pending

The complete scPerturb Adamson UPR source now passes native annotation and
streamed pseudobulk reconstruction. An independent HDF5/SciPy check verifies
all original datasets and every aggregate/QC count. No perturbation predictor
has been fitted or scored on this study. The [frozen protocol](PROTOCOL.md)
remains unchanged, SHA-256
`d6e73cda0945d0e8c9f86f2ea2ab30ebfe5e59f0d6edc90871a8b34e146fb297`.

## Source and identities

[Adamson et al., Cell 2016](https://doi.org/10.1016/j.cell.2016.11.048),
PMID 27984733, [GEO GSM2406681](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2406681),
UPR CRISPRi screen 10X010 in K562. Input:
`AdamsonWeissman2016_GSM2406681_10X010.h5ad`, 471,286,951 bytes,
SHA-256 `e70fcd49808cab8d724de8d5a332940911206e1c8ef44cc7b568d048ed795c85`.
The MD5 `2fa44ea61a8dd35742af618638ec65fc` matches the indexed primary
[scPerturb Zenodo 7041849 record](https://zenodo.org/records/7041849).
Direct Zenodo requests returned gateway errors; acquisition used the
[official scverse mirror](https://exampledata.scverse.org/pertpy/adamson_2016_upr_perturb_seq.h5ad)
identified by the [pinned pertpy loader](https://github.com/scverse/pertpy/blob/5bf4acaa45af16b49d248e4986fe394344fafdf5/src/pertpy/data/_datasets.py).
Direct current-release Zenodo metadata was not retrieved; no claim of release
1.4 identity is made. Download headers, hashes and failure logs are retained.

| Audited quantity | Complete source |
| --- | ---: |
| Cells | 65,337 |
| Genes | 32,738 |
| Stored/canonical nonzeros | 237,812,947 |
| Total integer UMI counts | 1,039,857,798 |
| Largest stored count | 3,449 |
| Explicit zeros / zero-library cells | 0 / 0 |
| Source guide categories | 114 |
| Cells with missing guide codes | 2,613 |
| Aggregate groups, including missing | 115 |

The source's `nperts` values are 0 for 2,613 cells, 1 for 101 cells and 2 for
62,623 cells. They are **not target multiplicities**: the pinned
[scPerturb preprocessing helper](https://github.com/sanderlab/scPerturb/blob/b69f72a070a92bcbaf41e7f9897b11598109ab48/utils.py)
derives this field by splitting the label on underscores. Guide names contain
such separators. The 101-cell `*` category also remains explicit.

`62(mod)_pBA581` (2 cells), `63(mod)_pBA580` (6,010) and
`Gal4-4(mod)_pBA582` (1,283) have not yet been verified against primary
control/guide records. They remain separate source-guide aggregates. No
negative-control assignment is inferred from the labels. Author demo code found
so far describes another experiment, and the public supplementary API reported
the requested article unavailable for that route. These missing identities
prevent the frozen response benchmark from proceeding to fitting.

## Native change and verification

Previously, annotation rejected this file at its 64 MiB whole-file read limit.
It now snapshots/hashes and publishes with 1 MiB buffers, admitting up to 1 GiB
sources and 2 GiB results. Existing axis/payload/storage checks remain. The
publication helper retains descriptor-rooted atomic no-overwrite behavior.

The source-bound annotation adds only `obs/numivivo_source_guide`, preserving
every existing label and replacing missing codes with a clearly named category
in that new column. The original `obs/perturbation` is untouched. No source cell
or gene is dropped. Native aggregation uses these source-guide identities and
an explicit pooled-replication-unresolved identifier, with no invented donors,
DE contrasts, control inference or gene-prefix pooling.

The independent reference multiplies an integer sparse group-membership matrix
by CSC blocks of 64 genes, auditing all entries. Only the 115-by-32,738 aggregate
is dense; the cell-by-gene array is never materialized.

On the physical M4 Pro (24 GiB, macOS 26.6 build 25G72, Swift 6.3.3):

- Release product build passed; frozen executable SHA-256
  `93986de58266a8cbc6f43f51537907ad1e39ba8a09a0bdbbfdfcdf95f9b72469`.
- All 30 original HDF5 datasets retain exact values, dtypes, shapes, attributes
  and compression/chunk settings. AnnData's backed reader accepts the output.
- All aggregate counts, per-cell total/detected/mitochondrial counts, feature
  identities and source-cell memberships match the independent reference exactly.
- Native pseudobulk verification reconstructs the complete report exactly.
- The existing annotation interoperability suite passed 26 checks; the new
  sparse-file admission fixture rejects 1 GiB + 1 byte without publishing output.
- Ten tests passed across rooted publication, interchange and target-kernel
  suites. The standalone focused publication test is included in those ten.
- Annotation took 2.53 seconds, 85,196,800 bytes peak RSS; native aggregation
  took 39.60 seconds. Test compilation overlapped aggregation, so these are run
  diagnostics, not a controlled CPU/Scanpy performance comparison.

Output H5AD: 475,168,719 bytes, SHA-256
`8ed5c635c9e5408bede2debc2104eb02ab7769ae4599eddd4296361ee4cbc13f`.
Native report SHA-256:
`79466c55cda1f4ffb47f689e836862010eeaab8d76f11819a4a5742196da23a6`.
The output receipt binds executable/runtime identity; a different build or OS
requires fresh receipts. HDF5 metadata allocations remain outside the fixed
copy buffer. Million-cell annotation, scientific prediction, independent donor
replication and Metal speed are not established by this check.

## Reproduce

Use the pinned original file and a built `numivivo` with
`NUMIVIVO_HDF5_LIBRARY` set to a compatible HDF5 library. Python reference tools
use NumPy, SciPy, h5py and AnnData; versions are in the evidence archive.

```sh
python prepare_ingestion.py --source original.h5ad --out prepared
numivivo singlecell-h5ad-annotate original.h5ad --plan prepared/annotation-plan.json --output annotated.h5ad > annotation-receipt.json
numivivo singlecell-h5ad-pseudobulk annotated.h5ad --plan prepared/pseudobulk-plan.json --output pseudobulk
numivivo singlecell-h5ad-pseudobulk-verify pseudobulk
python check_ingestion.py --source original.h5ad --annotated annotated.h5ad --annotation-receipt annotation-receipt.json --prepared prepared --pseudobulk pseudobulk --out check.json
```

The [evidence manifest](evidence/2026-09-10-native-ingestion/manifest.json)
hashes every retained file and records compressed logical-file hashes. Large
source/output H5ADs and executable remain external, with exact identities and
reconstruction inputs recorded; they are not embedded in Git. The source-guide
count reference and native report are retained. Next: verify primary guide and
control identities, group all guides for each target together, capture fixed GO
descriptors, then run and score the predeclared held-target protocol.
