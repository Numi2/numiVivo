# HIRISA complete-source experimental benchmark

The complete deposited GSE306664 release is acquired, audited and available as
one verified AnnData source: **1,612,594 unique cells, 18,082 genes,
3,845,991,249 nonzero count entries and 7,700,096,227 UMIs**. All 131 original H5
files are retained unchanged. Native execution, paired DE and held-out prediction
remain pending; this is source/interoperability evidence, not a completed
biological or native million-cell benchmark.

The [protocol](PROTOCOL.md) was frozen before count inspection, fitting or
prediction. Its SHA256 is
`24680681c348a92f75308225d6c162f20c7e9be8191c19a1339ab54c499ee10d`.
The source SOFT metadata and the first file's HDF5 hierarchy/dtypes were inspected
before freezing. No response signatures or DE outcomes selected this cohort.

## Deposited source and design

| Original enrichment | Libraries | Cells retained |
| --- | ---: | ---: |
| Bcell | 25 | 300,665 |
| Monocyte | 25 | 319,578 |
| NK | 25 | 289,717 |
| Tcell | 30 | 375,572 |
| Whole PBMC | 26 | 327,062 |

Five original donors contribute all populations. The metadata-only matcher
materializes sixteen enrichment/treatment contrasts. Fifteen have five donor
pairs; Bcell IFNg has four. GSM9205558 is a later-batch Bcell IFNg library with
no same-batch enriched control. Its cells remain in the source, with that donor
pair withheld from inference. Tcell controls from separate experiment batches
are matched separately. Whole-PBMC Fresh libraries are not used as 21-hour
unstimulated controls; culture_IFNa/culture_no_stim and technical pool identities
remain distinct for the later cross-preparation test.

These are deposited labeled 10x Flex counts. They include author DC and other
labels even though the paper describes removing such cells later. The paper's
curated analysis cohort and its reported UMAP population are not equated with
this complete release. Author labels and confidence scores are preserved as
annotations, not verified cell identities. Original cell_uuid is globally
unique across all files; no accession prefix hides duplicated cells.

## Exact interoperability checks

`acquire.py` verifies all 131 expected file sizes from the GEO file listing,
records SHA256 identities, retains source metadata, and freezes the control
matching without expression data. The original payload is 4,741,024,565 bytes.

`audit_sources.py` independently traverses each gene-by-cell CSC source in
1,024-cell blocks using h5py/SciPy, verifies sparse structure and integer counts,
and records exact per-library gene sums and per-cell UMI/detected-gene totals.
All 1,612,594 deposited n_umis and n_genes values agree. Global identity checking
uses a disk-backed SQLite table. The initial 26-file and 75-file partial audits
are retained; the final pass reuses 75 per-file count checks only after verifying
source, auditor and reference hashes, and reruns the global identity check.
Its 25.19-second final-pass time is therefore not a fresh full-count scan time.

`prepare_h5ad.py` retains all entries as CSR in their original cell/feature
order, with 64-bit offsets because the entry count exceeds signed 32-bit range.
It checks every copied count and index, preserves all 18 original observation
columns and adds seven explicit GEO design columns. AnnData string fields use
variable-length UTF-8 encoding. Axis-aligned feature annotations are in var;
remaining original feature tags/target-set structures remain in the unchanged
source H5 files. No cells by genes dense matrix is created.

The prepared H5AD is 6,117,413,997 bytes, SHA256
`0873e698ebf8770a54e6dba09724ffbeda5e1a67dbf24d4e223552d4aca7969c`.
Preparation took 217.12 seconds with 187,613,184 bytes maximum RSS on the local
Mac. This is Python preparation, not native NumiVivo throughput.

`verify_h5ad.py` then opens the entire result through AnnData 0.13.3.post0's
public backed reader and scans every sparse row. All per-cell UMI and detected-
gene totals and all 131 full gene-count aggregates match the independent source
references. It verifies original UUIDs and full axis sizes. The verification
took 95.42 seconds with 2,910,240,768 bytes maximum RSS; observation metadata is
resident despite the backed sparse matrix. This is an interoperability check,
not an end-to-end native/scverse performance comparison. The receipt's
allCellQCVerified flag refers to UMI totals and detected genes for every cell;
it does not assert other biological QC metrics.

## Native execution and next acceptance gates

The full source is transferring to the physical Mac mini for the unchanged
native executable at implementation 865f6f6a1963389f4fa4a475fe06b34712123c2c,
SHA256 `c4f7b08113215cb66413ab4ab62b868b57814d81332b71166f5def25a097c0ab`.
A pre-transfer capacity guard first rejected insufficient free space; that
rejection is retained. Verified unused caches and redundant transfer bundles
were cleared before retrying. The archived transfer snapshot is time-specific.
Inspect the live driver before retrying or asserting terminal status.

Code inspection shows the current native streamed source limits are one million
cell identities and one billion entries, below this complete input. The actual
full-source native baseline must finish before recording its precise failure.
Then repair the owning admission/storage paths, preserve failures, and qualify
complete publication/replay and every native aggregate/QC value. Metadata and
report residency remain unresolved even after sparse entry streaming.

Continue the frozen sixteen paired contrasts with native Wald/LRT/adjusted QL
and matched edgeR/limma/DESeq2 references, then all eligible held-out donor
prediction folds. Before fitting, freeze the existing native prediction options
and the separate PBMC annotation/transfer mapping. No production default,
biological calibration, million-cell full analysis, PCA/graph scale or Metal
speed claim follows from the current source checks.

## Reproduction and evidence

Use the pinned [requirements](requirements.txt). Each script uses the same root:

```sh
python acquire.py --root /absolute/path/hirisa
python audit_sources.py --root /absolute/path/hirisa
python prepare_h5ad.py --root /absolute/path/hirisa
python verify_h5ad.py --root /absolute/path/hirisa
```

Keep failed partial downloads/preparations for diagnosis; the scripts refuse
silent overwrites. `--available-only` makes partial source audits explicit and
does not qualify the complete release. Check [evidence/2026-09-10](evidence/2026-09-10)
with `python verify_archive.py evidence/2026-09-10`. The archive retains all
source hashes/URLs, design, per-file audits/references, command logs and receipts.
Original H5 files, the complete H5AD and the rebuildable global UUID index remain
external. Native results must be archived separately once their live run is
terminal; the current archive records pending status rather than a pass.

Sources: [GEO GSE306664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306664),
[author resource](https://apps.allenimmunology.org/aifi/resources/ifn-response/),
[experimental methods](https://apps.allenimmunology.org/aifi/resources/ifn-response/methods/),
[analysis methods](https://apps.allenimmunology.org/aifi/resources/ifn-response/analysis/),
and [AnnData format specification](https://anndata.readthedocs.io/en/stable/fileformat-prose.html).
