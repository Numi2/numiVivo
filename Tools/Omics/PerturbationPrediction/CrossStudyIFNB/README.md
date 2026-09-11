# Kang–HIRISA IFN-beta transfer

**Biological primary comparison: FAIL in both directions. Numerical execution:
PASS for all 26 folds and all 104 prediction vectors.** This experiment tests the
frozen donor-response baselines across two previously inspected studies. It adds
an actual cross-study prediction result, not untouched external validation.

All 2,651 Kang B cells and 119,513 HIRISA enriched-Bcell cells contribute to the
admitted cohorts: eight and five paired donors respectively. The fixed panel has
11,884 exact, unique source gene-symbol matches, out of 15,706 Kang and 18,082
HIRISA features. Every source feature remains in its library's normalization
denominator; absent genes are not padded with zero and aliases are not guessed.
Source Ensembl IDs, symbols, unmatched/ambiguous features and original group/cell
memberships are retained in the input mapping and cohort records.

## Complete results

Lower is better. These are equally weighted donor mean response RMSEs over the
entire shared panel in natural log1p(CPM) units. Each cross-study fit uses all
training-study donors. Matched within-study references exclude the query donor.
The [protocol](PROTOCOL.md) fixes the primary comparison before cross-study fitting:
ridge must beat both no-change and the cross-study training mean.

| Training → query | Query donors | No change | Cross mean | Cross median | Cross ridge | Within ridge | Primary |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| HIRISA → Kang | 8 | 1.220670 | 1.135572 | 1.135573 | 1.178674 | 1.132625 | FAIL |
| Kang → HIRISA | 5 | 0.341591 | 0.521068 | 0.453388 | 0.596829 | 0.107549 | FAIL |

Cross ridge is worse than cross mean for **all 13 query donors**. Its direction
mean is 3.80% worse for Kang and 14.54% worse for HIRISA. In Kang it improves on
no-change for all eight donors; in HIRISA it is worse than no-change for all five.
The simple HIRISA mean response has some transfer signal toward Kang, but neither
direction qualifies the context-ridge primary comparison. No model is tuned or
promoted from these scores. Full donor/method errors and correlations are archived.

## What this comparison means

[Kang's author manuscript](https://stacks.cdc.gov/view/cdc/79371/cdc_79371_DS1.pdf)
describes six-hour IFN-beta stimulation of PBMCs from lupus donors.
[HIRISA's deposited experiment](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306664)
uses healthy-donor enriched B cells, 100 units/mL IFN-beta for 21 hours, and fixed-RNA
probe measurements. Health status, duration, preparation, chemistry and study vary
together. Gene matching does not make these assays quantitatively interchangeable.
The result measures combined context transfer, not a causal effect of one difference.

Kang author-labeled B cells and the complete HIRISA Bcell-enriched libraries are
not asserted to be identical pure cell populations. Both studies' outcomes were
already inspected in earlier within-study work. Their donor identifiers are
retained with explicit study namespaces; no genotype-based participant linkage
was performed. There are thirteen biological donors, not 26 independent studies
or 122,164 independent experimental replicates. No calibrated intervals, new
patient outcome prediction, or RNA-to-phenotype link is established.

## Native implementation and verification

`VivoPerturbationPlan.responseFeatureIDs` optionally declares an ordered panel for
both response outputs and context-feature selection. Normalization uses the full
input library before projection, in both training and prediction. Every panel
gene must be present in each query feature dictionary; extra measured genes
contribute to its denominator. The default remains the strict identical-feature-
universe contract. `impliedCPMSum` is a panel subtotal when the option is present.

The actual CLI now accepts fit plans up to the owner's existing 2 MiB bound.
The initial real fit failed before model construction because its 128 KiB CLI
loader rejected this panel; that attempt and executable remain retained. The
fixed run uses the unchanged frozen inputs. Two standalone test compilation
attempts lacked the Testing framework/plugin paths; `test.sh` records the working
physical-Xcode invocation. These are retained software failures, not failed
biological fits or reasons to change the prediction protocol.

All 26 native fits, model verifications, predictions and prediction verifications
completed, plus one exact repeated prediction: **105 successful commands**.
Five Swift tests pass, including full-denominator behavior, feature reordering,
missing/duplicate/empty panel rejection and legacy plan encoding. Eight additional
actual CLI regression commands pass their declared success/rejection outcomes.
The default 18,082-gene HIRISA fit has exactly the historical numerical model and
all four prediction vectors; only transport identities differ.

Independent NumPy 2.5.3 / SciPy 1.18.1 / scikit-learn 1.9.0 reconstruction verifies
all native training count vectors, feature selection, moments, coefficients and
104 prediction vectors. Maximum prediction difference is **1.262e-13**;
maximum dual-coefficient difference is **4.285e-13**. Two complete scoring runs
produce identical checks, comparisons, scores and summaries. Native provenance
checks and numerical agreement do not turn the failed biological comparison into
a pass. This qualifies the scoped CPU Omics owner/CLI, not the entire application
or GPU prediction performance.

Physical execution: M4 Pro Mac mini, Swift 6.3.3, macOS 26.6. Final native executable
SHA256: `9fd1a3934ccfb921a47b08bf8728683641abe03ae8f012af968794e75a3381cc`.

- Input freeze: `3fde88150541b518f6bbd0e69c87f6d7fa35bafb7877b23703f9163cbbce8c5f`.
- Prediction freeze: `02ff87b65843a80d513caaf68cf7113f54439892ce60cd1764a47c5d1b5c2c27`.
- Scorer: `71ba6c8d58a295ab9c4a0ea1619f8697a4cfe3391609d68e1b291119b16241f2`.

## Reproduce

The transport consists of complete donor pseudobulk rows, explicitly not
individual-cell observations. `prepare.py` checks the prior native Kang report
identity and HIRISA's complete native/source count linkage, then compares every
selected HIRISA count vector with that linked source aggregate. The previously
qualified source studies remain the owner of original raw-cell ingestion.

```sh
python prepare.py --kang /path/to/kang-gamma-product \
  --hirisa /path/to/hirisa-study --out /new/inputs
bash Tools/Omics/H5AD/build.sh /new/runtime --with-cli
bash Tools/Omics/PerturbationPrediction/CrossStudyIFNB/test.sh /new/runtime
NUMIVIVO_HDF5_LIBRARY=/path/to/libhdf5.dylib python run.py \
  --inputs /new/inputs --binary /new/runtime/numivivo-omics --out /new/native
python score.py --inputs /new/inputs --native /new/native --out /new/score
```

Use the scripts in this directory and the repository-relative build/test paths
from the checkout root as appropriate. Python preparation requires AnnData,
NumPy, pandas and SciPy; scoring also requires scikit-learn. No Python model
produces the native predictions. `regression.py` additionally accepts the retained
historical HIRISA batch fold `005` to verify default full-universe behavior.

The [compact evidence archive](evidence/2026-09-11) stores original metadata,
freezes, scores, source/build identities, logs and original source code. Large
H5AD, NPZ, model and prediction payloads remain in
`/Users/n/numivivo-cross-study-ifnb-20260911` on `macmini`; scientific JSON and frozen
inputs also reside in `/Users/home/numivivo-cross-study-ifnb-20260911`. The archive
lists their exact identities. Restore those external payloads before replaying
historical bundles; the compact archive alone is not the complete dataset.
