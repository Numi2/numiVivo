# Complete paired RNA/ATAC count benchmark

The benchmark uses the complete filtered [10x Genomics PBMC granulocyte-sorted
3k Multiome release](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0),
processed by Cell Ranger ARC 2.0.0 and published 2021-05-03 under CC BY 4.0.
The publisher describes paired ATAC and gene-expression libraries from isolated
nuclei of one healthy human donor. This supplies paired assay measurements,
without independent donor replication or a perturbation comparison.

## Units and identities

[Cell Ranger ARC 2.0 documents](https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/analysis/outputs/feature-barcode-matrices)
RNA matrix entries as UMIs and ATAC peak entries as **cut-site counts**. NumiVivo
now has an explicit `cutSiteCount` unit restricted to accessibility assays.
No factor-of-two conversion or relabeling as fragments is applied. Other
accessibility sources can still declare fragment/read counts when appropriate;
the source's measurement definition determines that choice.

The same source barcode identifies each nucleus in both spaces. All source
features are retained: 36,601 genes and 98,319 peaks across 2,711 nuclei. Peak
identities are checked against the separate source `features/interval` column,
with explicit assembly GRCh38 and zero-based half-open coordinates, including
non-primary contigs. The assay namespaces and units remain distinct.

## Reproduce

```sh
curl -fL 'https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5' -o original.h5
python check_multiome.py --binary /absolute/path/numivivo --source original.h5 --out qualification
```

The input has 38,844,318 bytes and SHA-256
`5fbff5a4d85e0df345f6502e966ec787a8a4c429fd6b88a8772c43fd915cf3ff`.
`pbmc3k-multiome-plan.json` fixes the sample, feature namespaces, measurement units
and evidence provenance. The checker rejects a different source fingerprint.
No cell or feature selection is performed.

The independent oracle constructs a SciPy sparse matrix directly from the HDF5
source and checks Scanpy's complete-feature reader against it. It compares native
counts, identities, per-cell/per-feature totals and peak intervals exactly.
A second independent input constructs both AnnData modalities directly from the
original source arrays, writes CSC MuData, and imports that through the native
H5MU path. The resulting dataset and native H5MU must equal the 10x path byte for
byte. Repeated 10x imports must also produce identical bundles. Verification
reconstructs outputs from retained source and plan.

## Scope and resource bounds

The full source contains 24,511,186 nonzeros, exceeding the earlier 20-million
resident allowance. The common resident limit is now 32 million entries; source
files remain capped at 1 GiB and canonical dataset JSON at 512 MiB. The checker
records each command's elapsed time and maximum resident bytes through macOS
`time -l`. These measurements describe this pipeline and host, not a comparative
performance claim or million-cell qualification.

This benchmark qualifies measured paired count, identity and interval interchange.
It does not validate peak calling, fragment-file processing, TSS enrichment,
cell-type labels, peak-to-gene associations, RNA/ATAC integration, regulatory
prediction or biological mechanism. Those require their own methods and evidence.

## Qualified result on 2026-09-09

The release build and all 62 selected Swift tests in 17 suites passed with native
HDF5 enabled. All ten real-data product commands passed their expected outcomes,
including five deliberate rejections for incomplete modality mappings, misuse of
cut-site units on RNA, unknown units, assembly mismatch and excessive source entries.

| Assay | Features | Nonzero entries | Total measurement |
| --- | ---: | ---: | ---: |
| RNA | 36,601 | 5,218,473 | 11,786,194 UMIs |
| ATAC | 98,319 | 19,292,713 | 48,245,242 cut sites |

All counts, feature identities, paired barcodes and peak intervals matched the
source. The independent MuData route reproduced native dataset and H5MU bytes
exactly. Repeated native 10x imports reproduced the complete bundle exactly.

Peak resident memory across commands was 4,985,749,504 bytes (4.99 GB). The first
10x import used 3,291,955,200 bytes; H5MU verification had the largest footprint.
The measurements were taken on the laptop's Apple M4, macOS 26.6, with the release
executable built on the Mac mini. HDF5 was dynamically loaded from an isolated
HDF5 2.2.0 runtime. No Metal acceleration was used for these imports.

Evidence is retained in `evidence/2026-09-09-multiome`. The previous executable's
rejection of the correct `cutSiteCount` plan is retained, along with final test
and build logs, commands, fingerprints and real-data artifacts. Scanpy warns
about duplicate gene display names; stable feature IDs remain unique and no
renaming occurs. MuData's native-producer warning is retained as well.

To restore the native bundle, copy `bundle/original.h5`, `bundle/plan.json` and
`bundle/receipt.json`, then decompress its dataset JSON and H5MU files. To restore
the independently sourced H5MU bundle, decompress `independent.h5mu.gz` as
`original.h5`, use `h5mu/plan.json` and `h5mu/receipt.json`, and reuse the identical
native dataset/projection. Verification requires the recorded executable; a
new build must generate its own implementation-bound receipt.
