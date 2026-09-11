# Replogle 2020 UPR: complete inputs and confident cohort qualified

The complete deposited five-gemgroup experiment now passes native assay partition, RNA aggregation, reconstruction and independent AnnData/SciPy checks. **No predictor has been fitted or scored on this study.** The [frozen protocol](PROTOCOL.md) retains the existing fixed target-kernel method and all five technical conditions. This adds a separately collected experiment to the validation route while Adamson's unresolved identities remain open.

## Original information and assay separation

Source: [Replogle et al., Nature Biotechnology 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7416462/), DOI 10.1038/s41587-020-0470-y, PMID 32231336, [GSE146194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146194), deposited exp1–5 files under GSM4367979. The paper's Pilot UPR methods explicitly name `sgNegCtrl2` and `sgNegCtrl3` as non-targeting controls. This authority does not establish roles for the different Adamson construct IDs.

| Quantity | Complete source | Confident RNA cohort |
| --- | ---: | ---: |
| Cells | 40,997 | 32,829 |
| RNA features | 33,694 | 33,694 |
| Guide-capture features | 64, separately retained | Excluded from RNA denominators |
| RNA nonzeros | 129,837,664 | 105,512,809 |
| Guide-capture nonzeros | 1,913 | Separate assay |
| RNA groups | 165, including unassigned | 160 = 32 guides × 5 gemgroups |

All **129,839,577 original entries** and both complete feature axes are preserved. Original RNA UMI totals are 489,759,175; guide-capture totals are 2,483. Guide counts must not contribute to RNA library totals. Native projection produces separate RNA and guide H5ADs with identical complete observation axes; this is not full large-cohort H5MU/common-model qualification.

The original feature indices are unsorted within cells. Conversion preserves their order and verifies coordinate uniqueness within and across blocks. The original gene-by-cell Matrix Market coordinates are transposed into cell-by-feature CSR without a dense cell matrix. A semantic hash of every original coordinate/value triple is reconstructed from the converted H5AD. Independent checks compare every native partition coordinate/value and all axes.

The source assignment CSV has 34,406 rows; three barcodes are absent from the filtered matrix and are retained in an explicit outside-matrix record. Full barcode suffixes are used. Among source cells, 6,594 have no assignment, 702 fail guide coverage and 872 have multiple guide calls. The 32,829 included cells meet exactly `good_coverage=True` and integer-valued `number_of_cells=1`; both `1` and `1.0` CSV forms are validated without changing their meaning. These are guide-quality rules, not expression-based filters.

All five gemgroups and all 32 guide labels remain separate. They do not represent independent biological donors. Technical platform names, gene-level descriptor identities and missing/ambiguous mappings still require their declared authority before those labels enter prediction results. No guide is discarded to improve a score.

## Native change and validation

The previous projection owner rejected the real source at a fixed 100-million-entry cap. Sparse admission now uses the existing, charged element-visit allowance: default 500 million, maximum two billion. Source and transfer buffers remain bounded. Source snapshots retain the existing streamed owner's 64 GiB input bound; this experiment does not qualify files at that maximum.

A new optional `maximumOutputBytes` plan field permits positive output allowances up to 8 GiB. Omitting it retains the historical 1 GiB default and encoded plan. This experiment explicitly requests 2 GiB for the RNA projection; the guide projection uses the default. The allowance is private to each projection's HDF5 instance and governs allocation reservations and actual file-size checks. Full RNA requests with default storage or insufficient work reject without publication. No source subsampling or precision change was used to clear the former limit.

Validation on the physical M4 Pro, Swift 6.3.3 and HDF5 2.2.0:

- Both complete native assay partitions pass exact reconstruction.
- Every partition coordinate/value and all observation/feature identities match the original information.
- All 165 full RNA aggregate rows, cell QC values and memberships match the independent original-MEX reference.
- AnnData backed readers accept both complete partitions. SciPy reproduces all full RNA group sums and QC exactly.
- All 160 confident-cohort groups, 105,512,809 RNA nonzeros, memberships and QC match the independent SciPy reference; native replay passes.
- The projection regression suite passes 16 positive cases and 15 rejection cases. An explicit allowance preserves default H5AD output bytes; omitted optional fields remain absent from historical canonical plans.
- Two additional real-source checks preserve default output and insufficient-work rejection.

The executable compiles the actual single-cell product router and owners with `Tools/Omics/H5AD/build.sh <output> --with-cli`; it is a scoped build, not full-application qualification. Python/AnnData/SciPy are conversion and independent reference tools, not the native aggregation implementation. Build warnings and initial failed admissions/checker assumptions remain archived.

## Provenance and reproduction

Original compressed Matrix Market: **485,196,201 bytes**, SHA-256 `0684ccfa61d8460ecd9765afe020441de152e4ae2b91b9a5c6000221a9b1ea1a`.
Converted mixed-source H5AD SHA-256: `169b2a744e9188a05c4f5d3774bc43560c8b9316655673daef27eef3734b4acf`.
Full RNA report SHA-256: `3c07da62340f7f01f8b51022f4367785190dd803fb5214ed9363622c90d4461c`.

The [evidence archive](evidence/2026-09-11/manifest.json) retains source metadata, protocol freeze, author-code identity, plans, receipts, verification summaries and failure logs. Large source/derived matrices, aggregate reports, references and executables remain external under exact byte/hash identities. Do not replace those external artifacts with similarly named files or reinterpret historical receipts under another executable.

`prepare.py` converts and audits original GEO files. `run_native.py` partitions both assays and aggregates full RNA. `check_native.py` checks all source and native records. `check_scverse.py` performs the complete AnnData/SciPy check and freezes the confident-cell plan. `run_selected.py` executes and verifies that whole cohort. `check_large_admission.py` exercises real-source resource rejection. The study directory and exact binaries are recorded in evidence; destinations must be new.

Interactive GEO access and the PMC supplement download exposed challenges; neither was bypassed. Official public bulk GEO metadata/count downloads succeeded. The apparent spreadsheet response is retained as HTML, not accepted as a table. The original prediction-protocol hash preceded matrix acquisition; its source-declared RNA/guide distinction is recorded as an explicit pre-fit clarification.

Next: resolve exact descriptor identities and any required technical-platform naming, freeze every target fold and prediction, then run the separate scorer against all simple/shuffled baselines. Current GO knowledge may incorporate this and Adamson's studies. Success on an ingestion or replay check adds no biological prediction result, prospective target-selection evidence, new tissue context, calibrated uncertainty or clinical qualification.
