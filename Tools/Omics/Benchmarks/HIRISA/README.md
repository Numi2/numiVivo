# HIRISA complete-source experimental benchmark

The complete deposited GSE306664 release is acquired, audited and available as
one verified AnnData source: **1,612,594 unique cells, 18,082 genes,
3,845,991,249 nonzero count entries and 7,700,096,227 UMIs**. All 131 original H5
files are retained unchanged. All 48 native and 48 reference paired DE runs
are complete under the frozen execution specification. Complete native
ingestion/replay and independent cell/QC/aggregate and conditional DE numerical
checks pass. All 79 frozen donor-held-out prediction folds also pass native
replay and independent numerical checks, with [complete empirical results](PREDICTION_RESULTS.md).
This is not a completed biological or full million-cell analysis benchmark.

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

The full source reached the physical Mac mini. The unchanged baseline executable
SHA256 `c4f7b08113215cb66413ab4ab62b868b57814d81332b71166f5def25a097c0ab`
exited 65 with `HDF5 array exceeds admission limit`, in 3.152 seconds with
96,190,464 bytes maximum RSS. The preceding transfer-capacity rejection and this
actual native failure are both retained externally. The initial source archive
predates these results and still records a pending transfer.

The current reader admits two million source cells and four billion stored
entries, while retaining the separate aggregate and report bounds. Validated
sorted sparse segments bypass dictionary sorting; unsorted/duplicate segments
retain checked integer canonicalization. APFS snapshots now use independent
copy-on-write clones and hash the cloned bytes, falling back to descriptor copy
where cloning is unavailable. All 30 reader checks and four snapshot tests pass,
including independent source/snapshot writes. These are software gates, not
million-cell throughput qualification.

Complete ingestion succeeded against immutable debug executable SHA256
`d0cff005cc279fa7094232f1794bd1070b1c1d0c6d622c1ea3258d081b5e0a3e`,
in 1,630.69 seconds with 4,484,775,936 bytes maximum RSS. Its report is
526,518,683 bytes, close to the existing 512 MiB ceiling. Debug replay was
observed active at the original archive capture and subsequently passed in
1,667.41 seconds. The terminal receipt is in the inference archive; the earlier
archive's live observation remains unchanged. Metadata/report residency, PCA/graph scale and full out-of-core
qualification remain unresolved.

The same Omics source also built in release mode (271.29 seconds). Immutable
release executable SHA256
`0f676403a379783cbd6e48d75f9aa6e102090b4ca0b1eac7affe0270b6dea46e`
passed all 30 reader checks, then complete publication and reconstruction in
280.00 and 285.23 seconds, with 4,433,625,088 and 4,962,074,624 bytes maximum RSS.
Direct blockwise comparison proves the debug and release reports and plans are
byte-identical. Timings from these concurrently active runs are observations,
not a controlled native/scverse speed comparison. Both retained native replay
runs are now terminal successes.

`prepare_native_qc.py` independently reconstructs mitochondrial counts from the
complete backed sparse H5AD and rechecks every cell's UMI/detected-gene totals.
All 1,612,594 author mitochondrial-count annotations agree for the eleven mapped
mitochondrial genes. `verify_native.py` passed the independent native-result
acceptance gate: every report cell/QC object, all 131 source-library gene
aggregates, and the exact counts and original cell memberships in all sixteen
materialized inference cohorts agree. All 48 request hashes are bound to those
inputs. The audit took 421.35 seconds with 720,093,184 bytes maximum RSS on the
local Mac. It is separate from the successful release-native reconstruction.
The [native ingestion archive](evidence/2026-09-10-native) retains 490 members,
112,980,374 stored bytes, including the full compressed native report, exact
cohort inputs, cell/QC references, failed baseline, reader/snapshot gates and
frozen prediction folds. Its manifest SHA256 is
`43adc2bed11106286dea3b29f3779989204462e25d71cc30acc3df043d869a8b`.

The [inference execution specification](INFERENCE_EXECUTION.md) froze all 48
requests before fitting. Native Wald/LRT/adjusted QL use exact Python-prepared
source aggregates; the complete independent check now links those counts and
cell memberships to native ingestion. The matched edgeR/limma/DESeq2 runs and
native Wald/LRT/QL runs are all complete. Independent checks pass for all 48
native cases: counts/designs, offsets, fitted means, conditional scores,
likelihoods, information, effects, tails and BH arithmetic. QL tail checks use
the reported denominator degrees of freedom; they do not independently
requalify the moment/DF estimator on this experiment.

Across 144 native/reference comparisons, jointly tested effect-rank Spearman
correlation ranges from 0.954757 to 0.999573. Significant-call counts differ.
For example, the first Bcell IFNa contrast has 15,608 genes tested by all six
methods: BH 0.05 calls are 8,163/8,157/8,178 for native Wald/LRT/QL and
8,487/8,481/7,248 for edgeR QL/limma voom/DESeq2. Preserve rank-deficient support
statuses and each method's full tested family; do not interpret concordance as
ground truth, FDR control, power or interval coverage.

The [complete inference archive](evidence/2026-09-10-inference) contains 573
members and 750,221,134 stored bytes, including every native output and gene-level
numerical check, exact inputs, all reference tables/default DESeq2 filtering,
warnings, versions and the terminal debug replay. Its manifest SHA256 is
`c9d83cc81fd6ddd58ae598a06ede8c088297ec7e0be0447ddfeefe460b2cc2d1`.

The [prediction execution specification](PREDICTION_EXECUTION.md) and
`freeze_prediction.py` fix all 79 donor-held-out folds, the existing alpha-one
response model, training-only feature selection, and scoring against no-change
and training-response baselines. The [transport supplement](PREDICTION_TRANSPORT.md)
qualifies the native aggregate batch route without changing those model settings
or folds. All 79 donor-held-out fits completed; every output was frozen before
scoring. Complete raw-source reconstruction plus fitting took 271.56 seconds,
and native byte-exact replay took 276.32 seconds. Independent NumPy checks pass
for every fold, with maximum ridge prediction error 1.5954e-13. All four methods
are scored on the full and training-selected gene families. Mean, median and
ridge each improve contrast RMSE over no-change in 14 of 16 contrasts; ridge
improves on the simpler mean response in only four. Monocyte IFN-L1 and NK IFNg
remain worse than no-change for every learned baseline. See the full tables
and [prediction archive](evidence/2026-09-10-prediction), including the interrupted
provenance-correction attempt. No model settings were tuned from these outcomes.
The separate PBMC annotation/transfer mapping remains unfrozen. No production
default, biological calibration, million-cell full analysis or Metal speed claim
follows from these checks.

The 79-fold enriched-population experiment is complete at this stated scope.
Cross-preparation PBMC transfer and the million-cell PCA/graph/integration
pipeline remain separate, incomplete parts of the broader objective.

## Reproduction and evidence

Use the pinned [requirements](requirements.txt). Each script uses the same root:

```sh
python acquire.py --root /absolute/path/hirisa
python audit_sources.py --root /absolute/path/hirisa
python prepare_h5ad.py --root /absolute/path/hirisa
python verify_h5ad.py --root /absolute/path/hirisa
python prepare_inference.py --root /absolute/path/hirisa
python prepare_native_qc.py --root /absolute/path/hirisa
python freeze_prediction.py --root /absolute/path/hirisa
```

After complete native publication, copy its report/plan/receipt to the audit
host and run `python verify_native.py --root /absolute/path/hirisa --bundle
/absolute/path/native-full`. Preserve the separately executed native replay
status. For terminal inference cases, run `check_models.py --root
/absolute/path/hirisa/inference-inputs` and `summarize_inference.py --inputs
/absolute/path/hirisa/inference-inputs`. Their `--available` mode creates an
explicit partial checkpoint; it does not establish that all cases completed.
Pass `--native-verification /absolute/path/hirisa/native-independent-verification.json`
to the inference summarizer to include the separately checked native input link.

Keep failed partial downloads/preparations for diagnosis; the scripts refuse
silent overwrites. `--available-only` makes partial source audits explicit and
does not qualify the complete release. Check [evidence/2026-09-10](evidence/2026-09-10)
with `python verify_archive.py evidence/2026-09-10`. The archive retains all
source hashes/URLs, design, per-file audits/references, command logs and receipts.
Original H5 files, the complete H5AD and the rebuildable global UUID index remain
external. Verify the separate native ingestion archive with
`python verify_archive.py evidence/2026-09-10-native`. The original source
archive's pending-transfer record is historical; the newer archive records
completed release publication/replay and the independent full-source audit.
It retains the observed-live debug replay snapshot rather than inventing a
terminal result. Full DE output families and prediction remain separate gates.

The complete inference archive is verified with
`python verify_archive.py evidence/2026-09-10-inference`. Rebuild it from the
terminal external run using `archive_inference.py --root /absolute/path/hirisa
--output /absolute/path/new-inference-archive`. The new archive includes the
later terminal debug replay rather than changing the historical ingestion archive.

For prediction, preserve the original frozen fold manifest and transport receipt.
`prepare_prediction_batch.py --root /absolute/path/hirisa --repo /absolute/path/numivivo`
binds those folds to native group indices. Run the two native batch commands in
the [product documentation](../../PerturbationPrediction/NATIVE.md), then run
`score_prediction_batch.py --root /absolute/path/hirisa` only after all outputs
are frozen. `report_prediction.py` renders both complete contrast tables.
Archive with `archive_prediction.py` and verify every member with
`python verify_archive.py evidence/2026-09-10-prediction`. Large original source
and report payloads are shared through the earlier ingestion archive and external
source identity rather than duplicated per fold.

The [complete PCA qualification](PCA_RESULTS.md) retains all 1,612,594 cells and
compares every score to a frozen file-backed sparse reference. Native publication,
independent comparison and source reconstruction pass with unchanged numerical
settings. The selected cache is 2,263,542,272 bytes; QC JSON is 345,933,053 bytes.
These measurements drove shared storage-limit repairs. Metadata and fitted arrays
remain resident; complete graph results are recorded below and integration
remains a separate open gate.
Verify the full score/metadata/QC archive with `python verify_archive.py
evidence/2026-09-10-pca`; its chunk manifest checks complete reconstructed hashes
as well as every stored member. Earlier failed and interrupted runs are retained.

The [complete graph qualification](GRAPH_RESULTS.md) now passes native publication
and reconstruction on all 1,612,594 cells. Every returned distance and all
34,707,084 fuzzy edge records pass independent numerical checks. The fixed
2,048-query panel achieves mean strict recall 0.9997209821 against the complete
candidate universe; fifth-percentile recall is 1.0. Publication uses 6,963,750,537
metric evaluations and peaks at 4,360,978,432 resident bytes. The
[frozen protocol](GRAPH_EXECUTION.md) and 691-member archive retain source/runtime
identities, regressions, independent checks and transfer adjustments. Verify with
`python verify_archive.py evidence/2026-09-10-graph`. Clustering, integration,
biological preservation and Metal performance have their own acceptance gates.

Sources: [GEO GSE306664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306664),
[author resource](https://apps.allenimmunology.org/aifi/resources/ifn-response/),
[experimental methods](https://apps.allenimmunology.org/aifi/resources/ifn-response/methods/),
[analysis methods](https://apps.allenimmunology.org/aifi/resources/ifn-response/analysis/),
and [AnnData format specification](https://anndata.readthedocs.io/en/stable/fileformat-prose.html).
