# Restore original cell identities before evaluating Adamson predictions

The deposited scPerturb file has **3,634 incorrect or missing guide assignments**
relative to the original GEO cell records. Its preprocessing removes barcode
suffixes before joining metadata, then drops duplicated metadata indexes.
The 65,337 original full barcodes are unique; their bare prefixes have 2,539
duplicates across ten numeric suffix groups. AnnData-generated suffixes on the
deposited index do not recover those original identities.

This diagnosis is reproduced from the [pinned preprocessing code](https://github.com/sanderlab/scPerturb/blob/b69f72a070a92bcbaf41e7f9897b11598109ab48/dataset_processing/scripts/AdamsonWeissman2016.py),
the deposited file, and original public GEO records. We reconstruct both the
renamed observation index and the incorrect first-duplicate guide join exactly.
The GEO barcode order maps to all deposited matrix rows in that reconstruction.
No model scores or expression effects were used to detect or repair the join.

| Compared with original GEO records | Cells |
| --- | ---: |
| Labels recovered from missing deposited assignments | 2,535 |
| Existing labels replaced by a different original guide label | 1,097 |
| Deposited labels with no original assignment | 2 |
| Total changed assignments | 3,634 |
| Missing guide records after restoration | 80 |

The affected cells contain **55,413,685 UMI counts**. The historical native
ingestion check correctly reproduced the deposited data, including its erroneous
labels. It remains evidence of file/count interoperability; its guide aggregates
must not be used as authoritative experimental assignments for prediction.

## Original records and native restoration

Retrieved through GEO's public data endpoint:

- [GSM2406681 barcodes](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2406nnn/GSM2406681/suppl/GSM2406681_10X010_barcodes.tsv.gz),
  249,427 bytes, SHA-256
  `1e0d820343d0c6e17bdab4fef96f4446e223b94861942ccf3d7225818e009836`.
- [GSM2406681 guide records](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2406nnn/GSM2406681/suppl/GSM2406681_10X010_cell_identities.csv.gz),
  808,896 bytes, SHA-256
  `8b40be7a2280c1713bf5a1eb828aad46e70ad3a7000b5f5a1322e51e06c4cf7f`.

`restore_geo_identities.py` checks both hashes and the original H5AD hash,
reconstructs the source preprocessing, and supplies a source-bound native
annotation plan. Native Swift/HDF5 adds nine `obs/numivivo_geo_*` columns:
original barcode, guide label, guide/suffix sample identity, suffix group,
good-coverage flag, number-of-cells field, guide read count, guide UMI count and
guide coverage. It also appends the native provenance journal. All original
fields and X remain unchanged.

The corrected mapping explicitly uses the restored barcode and sample columns.
There are 1,106 guide/suffix sample combinations, pooled into 115 guide groups.
Numeric barcode suffixes are retained as technical GEM-group identifiers, with
no assertion of independent donor or biological replication. No cells are
filtered. Source flags include 59,219 good-coverage records, 6,038 false records
and 80 missing records; these remain available for a later explicit cohort rule.

Preserving this metadata exposed a CLI/engine inconsistency: the 193,098-byte
pseudobulk mapping failed the CLI's 128 KiB generic document reader although the
native engine admits 2 MiB plans. The pseudobulk CLI now reads up to 2 MiB, matching
publication and verification. The initial rejection is retained. A 2 MiB + 1
byte plan still fails before reading the matrix and publishes no output.

## Verified result

- The release product rebuilt on the physical M4 Pro. New executable SHA-256:
  `a4b44ad3e191e29bc80b25d0723002280064f108b08d14be65c1625c709615b4`.
- Native annotation used the previously qualified executable
  `93986de58266a8cbc6f43f51537907ad1e39ba8a09a0bdbbfdfcdf95f9b72469`.
  The two steps retain their actual, distinct implementation receipts; no
  receipt has been relabeled after the CLI rebuild.
- All 30 original HDF5 datasets retain exact values, dtypes, shapes, attributes
  and storage settings. AnnData opens the restored file in backed mode.
- Every restored metadata value matches the original GEO files, including
  nullable quality fields. All source cells and feature identities remain.
- Native aggregate counts match independent integer SciPy aggregation exactly.
  Per-cell quality and global per-gene totals match the previous source audit.
  Every original barcode, sample and technical-batch membership is verified.
- Native pseudobulk verification reconstructs the complete new report exactly.

Corrected H5AD SHA-256:
`a33bab097da97dc46414473d1761675940c59df8a53ee8517f60f98f950c334e`.
Corrected report SHA-256:
`9fd456f892ad2fc037194e154da640ea6293bbae79427ac75017ad030332ae4a`.
Counts remain 65,337 cells, 32,738 features, 237,812,947 nonzeros and
1,039,857,798 UMIs. Annotation took 3.32 seconds with 161,251,328 bytes peak RSS;
aggregation took 41.21 seconds with 544,473,088 bytes peak RSS. These describe
this bounded run, not million-cell, biological or GPU qualification.

## Frozen GO descriptors and remaining biological decisions

All 90 gene-like source-guide prefixes match unique original Ensembl features;
those identities and candidate names are unchanged by the GEO restoration.
`prepare_descriptors.py` captured current MyGene build `20260906` without reading
expression outcomes. Independent reconstruction verifies 2,932 retained
annotation records, 895 distinct direct GO terms and 88 supported candidates.
IER3IP1 has no usable terms; the source Ensembl ID for TIMM23 returns 404. No
current-symbol or alias substitution repairs missing identities. All supported
candidates share at least one term with another candidate. Pairwise set Jaccard
and independent binary-incidence calculations agree exactly; the minimum kernel
eigenvalue is 0.2091193544. Six corruption checks also reject modified, missing,
unlisted or rehashed-inconsistent evidence.

One retained MANF annotation, GO:0036500 (ATF6-mediated unfolded protein
response), uses IMP evidence and directly cites the Adamson paper, PMID 27984733.
This is explicit knowledge overlap. The frozen protocol uses current annotations
and retains that record; it does not establish temporal independence from the
study. No post-score annotation selection has occurred.

The original GEO CSV confirms the cell-to-guide labels but does not resolve the
negative-control constructs `62(mod)_pBA581`, `63(mod)_pBA580` and
`Gal4-4(mod)_pBA582` or the interpretation of all coverage fields. Controls remain
unassigned; no response model has been fitted or scored. Before fitting, resolve
those experimental definitions and apply the predeclared unambiguous-condition
rule, retaining exclusions. The [prediction protocol](PROTOCOL.md) is unchanged.

## Reproduce

```sh
python restore_geo_identities.py --source original.h5ad --barcodes GSM2406681_10X010_barcodes.tsv.gz --guides GSM2406681_10X010_cell_identities.csv.gz --out restored
numivivo singlecell-h5ad-annotate original.h5ad --plan restored/annotation-plan.json --output restored/annotated.h5ad > restored/annotation-receipt.json
numivivo singlecell-h5ad-pseudobulk restored/annotated.h5ad --plan restored/pseudobulk-plan.json --output restored/pseudobulk
numivivo singlecell-h5ad-pseudobulk-verify restored/pseudobulk
python check_geo_restoration.py --source original.h5ad --barcodes GSM2406681_10X010_barcodes.tsv.gz --guides GSM2406681_10X010_cell_identities.csv.gz --restored restored --old-prepared prepared --out check.json
```

`prepared` is the independent complete-source reference from the prior ingestion
step. Native commands require the HDF5 library; Python is the independent metadata
adapter/reference and is not called by native annotation or aggregation.
The [evidence manifest](evidence/2026-09-10-geo-restoration/manifest.json) preserves
raw GEO records, GO captures, plans, original and failed execution receipts,
corrected aggregate/reference counts and validation results. Large H5ADs and
executables remain external with exact hashes. Redundant task-owned H5AD copies
were removed only after hash/handle checks; both historical and corrected source
bundles remain retained on the Mac mini.
