# Full-source pancreas benchmark

The added dataset is the complete human portion of
[Baron 2016, GSE84133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133).
All 8,569 deposited human cells and all 20,125 genes are retained. The four
human matrices contain 16,171,764 nonzero integer entries and 49,942,062 UMIs,
spanning four donors and 14 author-assigned cell types. The two separate mouse
matrices are not combined with human genes or counted as human donors.

## Source and design

`prepare_baron.py` pins both the original GEO tar archive and family SOFT
metadata. It reads one CSV cell row at a time, validates nonnegative integer
counts, and verifies identical gene axes before constructing sparse CSR H5AD.
It never intersects genes, filters cells, rounds expression, or densifies the
cell-by-gene matrix. Python preparation/reference code still holds the full
sparse matrix in memory; it is not the native streaming implementation.

The source cell identifier becomes the AnnData observation identifier and native
barcode identity. The original sequencing barcode remains a separate
`source_barcode` observation column, retained in the unchanged H5AD snapshot.
The source `assigned_cluster` labels are retained verbatim as `cell_type`.
These author computational annotations are comparison references, not independent
experimental ground truth or authoritative automatic NumiVivo labels.

GEO explicitly identifies human1–3 as non-T2D and human4 as T2D. The mapping
preserves that distinction. It does not invent paired conditions, replicated
disease contrasts, or technical batches. `batchID` is `unreported`; donor identity
comes from the deposited human sample and its associated donor description.

Primary downloads:

- [Original count archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84133/suppl/GSE84133_RAW.tar)
- [GEO family metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84133/soft/GSE84133_family.soft.gz)

Full SHA-256 values, per-member hashes, source metadata fields, the prepared H5AD
hash, per-cell-type donor composition and qualification results are in
[`evidence/2026-09-09-pancreas`](evidence/2026-09-09-pancreas/).

## Verified native scope

| Check | Result |
|---|---|
| Every original CSV row vs prepared sparse H5AD | Exact values and identities |
| Every native cell total and detected-feature count | Exact source and Scanpy agreement |
| Donor/cell-type pseudobulks | All 56 exactly match sparse source sums |
| Source cell membership | Every cell occurs exactly once |
| Native bundle reconstruction | Passed |
| Repeated native publication | Identical receipt |
| Resident count import | Controlled size rejection, no output |
| NB disease contrast, beta cells | Controlled insufficient-replication rejection, no output |

The native runtime scans HDF5 arrays in bounded slices and keeps metadata, cell
QC and pseudobulk aggregates resident. It processed all 16.17 million nonzeros
without raising the five-million resident count limit. Independent verification
compares every original count row and all native aggregates, not only shape or
hashes. No new native source was required: the tested binary is the published
`6a8658c` build, with exact binary/report hashes in `source-state.json`.

```sh
python Tools/Omics/Benchmarks/prepare_baron.py --archive /data/GSE84133_RAW.tar --soft /data/GSE84133_family.soft.gz --out /new/baron-prepared
python Tools/Omics/Benchmarks/check_baron_cli.py --binary /path/to/numivivo --prepared /new/baron-prepared --out /new/baron-native
python Tools/Omics/Benchmarks/check_baron_stream.py --archive /data/GSE84133_RAW.tar --prepared /new/baron-prepared --bundle /new/baron-native/bundle --out /new/baron-comparison.json
```

The first harness expected the wrong generic substring for the resident-size
error. Native execution correctly rejected it with exit 65. The failed expectation
and final corrected harness evidence are retained.

## Muraro is not an integer-count source

The original `GSE85241_cellsystems_dataset_4donors_updated.csv.gz` contains
3,072 cells and 19,140 features. All 12,442,034 nonzero values are fractional;
examples include 1.0019582262109 and 6.07143081403291. The entire table was
audited. It cannot qualify the current integer-count path and is not rounded or
inversely transformed to fabricate UMI counts. Its continuous expression values
require a separately typed assay path. The Hemberg-hosted author annotation URL
returned HTTP 403, retained in the evidence. Cell-name prefixes are recorded only
as prefixes, not promoted to verified donor metadata.

```sh
python Tools/Omics/Benchmarks/audit_muraro.py --source /data/GSE85241_cellsystems_dataset_4donors_updated.csv.gz --out /new/muraro-audit.json
```

## Remaining implementation and qualification

Full-cohort normalization, HVG selection and PCA now pass through the
[streamed H5AD reduction owner](../Reduction/STREAMING.md), with all 8,569 cells
and 20,125 features contributing to QC, normalization and feature selection.
The 2,000-gene PCA agrees numerically with Scanpy and reconstructs from the
archived source. This extends the earlier count-only evidence; that historical
evidence directory retains its original qualification boundary.

Baron also exposes a design limitation: native donor integration currently
rejects the disconnected donor/condition design, and the single T2D donor cannot
establish a replicated disease effect. Relabeling all conditions as one condition
would hide this fact. A future cell-type benchmark may evaluate mapping from
three non-T2D donors to the held-out T2D donor, with training-only preprocessing
and explicit rare-type support, but that has not been executed here. Broader
integration must measure cell-type and within-type state preservation without
claiming disease preservation from an unidentifiable contrast. Streaming downstream graphs/integration,
multi-cell-type integration, learned reference mapping, perturbation prediction
and million-cell performance remain open.
