# HIRISA preparation-transfer execution

The [fixed protocol](CONTEXT_TRANSFER_PROTOCOL.md) now defines both enriched-to-PBMC
and PBMC-to-enriched IFNα response prediction, withholding the query donor from
all training preparations. Sixty cross-preparation folds and sixty matched
within-preparation references are frozen. **Predictions have not yet been fitted
or scored.** Native aggregation and replay are the next execution gate.

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

After all six native bundles pass, bind the frozen fold roles to native group
indices, run the unchanged native response owner, freeze all predictions, then
perform independent numerical reconstruction and score all outcomes. No native
count or prediction qualification is inferred from the reference preparation.
