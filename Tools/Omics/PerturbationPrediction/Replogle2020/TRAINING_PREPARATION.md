# Replogle UPR: native context preparation and count-fold freeze

**All five native training inputs and 150 guide-level count exclusions pass.
At this preparation stage, no model was fitted and no outcomes were scored.**
The subsequent [complete prediction experiment](RESULTS.md) is reported separately.
This completes control pooling and context isolation for the
[previously qualified complete cohort](README.md), under the unchanged
[prediction protocol](PROTOCOL.md). Experimental guide-to-gene descriptor
identities were unresolved at this stage and are now addressed by the separate
[identity audit](IDENTITIES.md). Count preparation cannot establish prediction
accuracy or substitute for that identity requirement.

## Complete cohort and controls

The native `singlecell-h5ad-pseudobulk` owner reads the retained RNA H5AD with
exactly the existing 32,829-cell confidence selection. It pools the raw counts
of the paper-defined `sgNegCtrl2` and `sgNegCtrl3` controls **within each original
gemgroup**. Their original sample IDs and cell memberships remain in the pooled
report, and the earlier 160-group audit retains each guide separately. Every
non-control guide remains its own condition; no expression-based filter is added.

| Original technical group | Confident cells | Pooled control cells | Pooled control RNA UMI | Complete selected RNA UMI |
| --- | ---: | ---: | ---: | ---: |
| gemgroup-1 | 8,006 | 1,763 | 22,740,829 | 99,075,431 |
| gemgroup-2 | 5,732 | 1,209 | 15,090,890 | 66,623,162 |
| gemgroup-3 | 5,736 | 1,528 | 18,884,523 | 66,943,446 |
| gemgroup-4 | 7,075 | 1,818 | 23,317,355 | 85,943,308 |
| gemgroup-5 | 6,280 | 1,085 | 14,071,274 | 75,455,095 |

The result has 155 rows: one pooled control plus 30 original guide conditions
per gemgroup. All 33,694 source RNA features and all 105,512,809 selected RNA
nonzeros remain accounted for. Technical names are not inferred from barcode
suffixes or expression. The original gemgroup labels are retained; neither
donor IDs nor biological replication are invented.

## Native preparation and independent checks

Conditions are explicitly named `gemgroup-N|guide` or
`gemgroup-N|pooled-control`. This uses the existing native mapping and selection
interfaces: each condition resolves to exactly one aggregate, and every selected
row shares its native biological-context fields. No Swift or prediction-method
change was needed.

The same actual single-cell executable used for cohort qualification executes
the new aggregate, its full reconstruction, and five
`singlecell-composition-prepare` commands. Each training file contains its
context's 31 rows, with the control first and original guide IDs sorted. Source
and aggregate-report fingerprints are checked. The driver independently verifies:

- Every gene count in all 155 native rows against the retained original-count
  SciPy reference, summing only the two declared controls where required.
- Every original source-observation index, selected cell QC, group membership,
  source sample ID, technical group and count-feature identity.
- Every value in all five native training inputs against those verified rows.
- Native rejection of a plan that takes a target from another gemgroup. It
  returns exit 65 with `composition conditions cross biological contexts`,
  and publishes no training file.

The reference uses condition-by-gene arrays and single-row comparisons; it never
allocates a dense cell-by-gene matrix. Existing input qualification is reused
under exact hashes rather than downloading or reconverting the cohort.

## What is frozen

`count-folds.json` retains every 30 × 5 guide exclusion. Each record names its
held guide, the other 29 training guides, the selected rows and the canonical
selected-training SHA-256. Changing every stored count in the withheld row to
one leaves the selected training object identical for all 150 exclusions.
The count-only subset contract matches the existing Norman native-fit driver.

These are **count-selection plans**, not fitted models or validated gene-level
folds. The retained context inputs contain all 30 guide rows for reconstruction;
a future fitter must receive only the recorded subset. Before fitting, verify
experimental gene identities, capture the exact GO descriptors and freeze the
descriptor/fitter/query inputs. Freeze every prediction before the separate
scorer reads held responses. Unresolved identities, unsupported annotations and
all five technical-group results must remain explicit. No endpoint or model
setting has changed, and no favorable fold has been selected.

## Provenance and reproduction

Executed on the physical M4 Pro against owner base
`72af0e8b02f57685605e077e180eb29df0866176`, using the retained Swift 6.3.3
scoped single-cell product and HDF5 2.2.0. This checks the existing native
preparation behavior, not a full application build or GPU prediction.

- Executable SHA-256: `36df30cd11afffb9bca0f88369af2dca3260759ddd9945c2c02d36073b1574db`.
- RNA H5AD SHA-256: `7170656afc5fd23fb7ee749467e5cb60832cc1766100aba8c99ec1ea96bfa17d`.
- Earlier separate-guide report SHA-256: `4c4142c29304d3f480bf5165232ad83d84d9acea0a8f536d604a9e4d62adec54`.

The [preparation archive](evidence/2026-09-11-training/manifest.json) retains the
driver, plans, all command logs, pooling audit, count-fold freeze and checks.
The check receipt records exact byte/hash identities of the external training
files and complete aggregate bundle metadata. The large H5AD source remains
under the qualified RNA identity above. Retained full reports and five training
files live in `macmini:/Users/n/numivivo-replogle2020-20260911/training-preparation`.

To reconstruct with the exact retained inputs and executable, choose a new output
directory with sufficient space for its source snapshot and outputs:

```sh
/Users/n/numivivo-replogle2020-20260911/anndata-env/bin/python \
  Tools/Omics/PerturbationPrediction/Replogle2020/prepare_training.py \
  /Users/n/numivivo-replogle2020-20260911 \
  --out /path/to/new-training-preparation
python3 Tools/Omics/Benchmarks/HIRISA/verify_archive.py \
  Tools/Omics/PerturbationPrediction/Replogle2020/evidence/2026-09-11-training
```

No downstream phenotype model, uncertainty coverage, independent laboratory,
prospective target selection or new tissue context is qualified by this result.
