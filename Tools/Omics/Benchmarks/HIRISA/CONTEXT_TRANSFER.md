# HIRISA preparation-transfer execution

The [fixed protocol](CONTEXT_TRANSFER_PROTOCOL.md) tests enriched-to-PBMC
and PBMC-to-enriched IFNα response prediction, withholding the query donor from
all training preparations. All **120 frozen folds completed**: sixty cross-preparation predictions and
sixty matched within-preparation references, spanning 705,365 selected cells and
all 18,082 source genes. All six native count bundles and prediction batches
pass independent numerical checks and native replay. Every output was frozen
before scoring, and two scoring executions produced identical result bytes.

All three learned methods beat no-change in **all twelve cross-preparation
contrast means**. Ridge meets the fixed
primary gate—lower all-gene RMSE than both no-change and cross-preparation mean—in
**3/12 contrasts**. Cross-preparation ridge is worse than matched within-preparation
ridge in all twelve contrasts, showing a consistent transfer penalty. No model
was selected or tuned from these results.

## Measured prediction results

Response RMSE is in natural-log(1+CPM) units; lower is better. Each entry is the
equally weighted mean of five held-out donors, over every source gene. Within
references use the same query and target, with the other donors from the query
preparation. Training-selected panels are secondary and differ between
preparations, so cross-minus-within comparisons use the all-gene family.

| Lineage / direction | No change | Cross mean | Cross median | Cross ridge | Within ridge | Ridge primary gate |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| B / enriched → PBMC | 0.297432 | 0.244057 | 0.244563 | 0.243124 | 0.196586 | pass |
| B / PBMC → enriched | 0.274059 | 0.186673 | 0.186943 | 0.189470 | 0.090472 | fail |
| Mono / enriched → PBMC | 0.461724 | 0.326023 | 0.334001 | 0.340340 | 0.245194 | fail |
| Mono / PBMC → enriched | 0.396165 | 0.310643 | 0.314633 | 0.296534 | 0.204975 | pass |
| NK / enriched → PBMC | 0.338269 | 0.253268 | 0.253980 | 0.253762 | 0.227802 | fail |
| NK / PBMC → enriched | 0.273513 | 0.183077 | 0.181050 | 0.188264 | 0.091108 | fail |
| CD4-T / enriched → PBMC | 0.282604 | 0.202605 | 0.203552 | 0.202752 | 0.156278 | fail |
| CD4-T / PBMC → enriched | 0.243838 | 0.170111 | 0.170263 | 0.174027 | 0.102762 | fail |
| CD8-T / enriched → PBMC | 0.434316 | 0.304924 | 0.306376 | 0.311911 | 0.248010 | fail |
| CD8-T / PBMC → enriched | 0.292021 | 0.256813 | 0.256802 | 0.261066 | 0.173676 | fail |
| other-T / enriched → PBMC | 0.422090 | 0.376468 | 0.379586 | 0.377552 | 0.357765 | fail |
| other-T / PBMC → enriched | 0.372447 | 0.333448 | 0.329929 | 0.331042 | 0.320424 | pass |

Across the sixty individual cross-preparation folds, mean, median and ridge are
worse than no-change in 4, 4 and 4 folds, respectively.
All fold scores, both feature families, implied CPM totals, clipped-gene counts
and twelve matched cross-minus-within comparisons remain in the
[complete execution archive](evidence/2026-09-11-context-transfer/manifest.json).
These are conditional expression point estimates; the annotation and
experimental-context limitations below still apply.

## Complete membership and reference counts

The frozen mapping uses exact author `celltype.l1` labels in the corresponding
enriched preparation and cultured PBMCs. The freezer directly verified the
original labels and UUIDs against all 131 deposited HDF5 files, and binds the
unchanged full source SHA-256
`0873e698ebf8770a54e6dba09724ffbeda5e1a67dbf24d4e223552d4aca7969c`.

| Author lineage | Cells selected across both preparations |
| --- | ---: |
| B | 132,026 |
| Mono | 157,533 |
| NK | 128,537 |
| CD4 T | 221,208 |
| CD8 T | 55,164 |
| other T | 10,897 |
| Total | 705,365 |

Every one of the original 1,612,594 cells has a frozen membership or exclusion
entry. The other 907,229 cells remain in the complete source and metadata ledger;
they are outside the fixed IFNα/cultured-control and annotation mapping. No count
threshold, DE outcome, annotation-confidence cutoff or model score selected cells.
The complete 18,082-gene feature universe is retained.

Nine original PBMC donor/batch/pool pairs contribute eighteen cultured libraries.
Raw technical-library counts are combined only within the same donor,
preparation, lineage and condition after pool matching. Enriched libraries use
their original matched control. Each lineage has twenty donor/preparation/
condition aggregates. All requested folds meet the ten-cell criterion; the
smallest aggregate has 178 cells. Pools never become additional donors.

The independent sparse reference verifies the frozen membership, all training
donor exclusions, accession separation and paired pool identities. It sums
**1,632,818,461 selected sparse entries and 3,285,105,717 UMIs** into 120 integer
aggregates, with a separate exact row-total check. Every library total is positive.
It processes bounded CSR blocks and materializes aggregate-by-gene arrays only.
The reference pass took 22.33 seconds with peak RSS 1,687,224,320 bytes on the
local host; this is Python preparation, not native or GPU performance evidence.

## Native and numerical qualification

The unchanged physical-M4-Pro release executable has SHA-256
`70a5bdcb258f35f947811bb7b5f2ddc777dfef19ce26880a74de818c8586ac95`.
All runtime source hashes were checked before execution. Each lineage's full
original-source aggregation, exact integer/membership comparison and native
replay passed. Prediction publish/replay then passed for all 120 requested folds.
All 120 excluded-count mutation checks produced the same projected inputs as
the native receipts. Independent reconstruction of fitted parameters and
expression vectors has maximum absolute difference **1.162e-13**. Implied CPM
totals were checked separately with the declared numerical tolerance.

The global prediction output freeze SHA-256 is
`17d298820cf2f4618bc9ab049468aa27baa95f46b2b5d1f3157018bfc5a6c2e8`.
The scored result SHA-256 is
`659181702f496192c678520e91157b7406a8a22a3197aaa94b843851ae7b4734`.
These checks establish provenance and numerical execution; empirical prediction
quality is the separate all-gene comparison above.

## Provenance and interpretation

The protocol was declared after the original enriched-population predictions
and integration diagnostics were known, before preparation-transfer fitting or
PBMC response scoring. It is a new held-out endpoint in an inspected study, not
an independent-study or prospective replication. Source counts may be prepared
for verification, but held-out treated rows are excluded from fit/query inputs
and predictions must be frozen before scoring.

Author labels are predicted annotations. Treated-cell labels define the scoring
strata and can depend on treatment. Preparation, experimental batch and culture
composition also differ together. These folds therefore test conditional average
RNA response across combined contexts; they cannot isolate a causal preparation
effect or establish prospective cell identity, unseen tissues or clinical utility.
All twelve lineage/direction contrasts, negative outcomes and unavailable folds
must be reported. A secondary feature panel cannot replace the all-gene endpoint.

Freeze SHA-256:
`408a2633f02a4c1586b274ce35fa536c2f93ecc739005944c0a51ffa17999b6a`.
Protocol SHA-256:
`53284b39b8d214d0321aabc3f01be9f715af33f92f5126733b6ed7fc54067434`.
The [archive manifest](evidence/2026-09-10-context-transfer-preparation/manifest.json)
retains all six complete selection plans, the compact original-row membership,
all folds and exclusions, aggregate identities and independent integer counts.

## Execute

```sh
python freeze_context_transfer.py --root STUDY --out STUDY/context-transfer
python prepare_context_reference.py --root STUDY --out STUDY/context-transfer/reference
python run_context_native.py --study STUDY --repo REPO --runtime FROZEN_RUNTIME --python REFERENCE_PYTHON --hdf5 HDF5_LIBRARY
python run_context_predictions.py --root STUDY/context-transfer --repo REPO --runtime FROZEN_RUNTIME --python REFERENCE_PYTHON --hdf5 HDF5_LIBRARY --count-coordinator EXISTING_PID
python score_context_predictions.py --root STUDY/context-transfer --out STUDY/context-transfer/scores.json
python archive_context_transfer.py --root STUDY/context-transfer --out ARCHIVE
python verify_archive.py ARCHIVE
```

The first two commands require the original local source files and Python with
NumPy, SciPy and h5py. The native driver uses the existing qualified full release
executable, checks all runtime source hashes and runs six complete lineage
aggregations sequentially. Each retains the full original source snapshot and
source-row selection. `check_context_native.py` compares every native count,
feature, group and original membership before native replay. The six plans
partition execution, not the scientific target; every selected cell is retained.
New output paths are required, and a terminal failure is preserved for repair.
Do not relaunch a live driver or restart completed phases after an observation
timeout.

The prediction driver waits for the existing count coordinator, verifies its
terminal receipt and maps each frozen fold to exact native group indices. It
runs at most two lineages concurrently, with a three-GiB free-space gate before
each phase. It requires a complete output freeze and native replays before
scoring. The count PID is only a liveness check while its terminal receipt is
absent; completed count work is reused without restarting.

For every fold, the preparer changes every stored count in excluded aggregates
in a temporary copy. The projected training and control hashes must remain
unchanged and then match the native fit/query receipts. The scorer independently
reconstructs feature selection, scaling, ridge fitting and all four predictions
from verified raw aggregates before comparing them with held-out outcomes.
The scorer also reproduces all 79 previously published HIRISA folds with exactly
identical metrics; that regression uses already known outcomes, not the new
preparation-transfer targets. All twelve contrasts require all five donor folds;
incomplete contrasts retain unavailable means rather than averaging successes.

The execution archive stores every model, prediction, native source report and
receipt in bounded compressed chunks. Repeated source files have explicit
restoration mappings; the unchanged full original H5AD remains externally
retained under its exact hash. The earlier preparation archive stays immutable.
