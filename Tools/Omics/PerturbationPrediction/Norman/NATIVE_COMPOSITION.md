# Native condition-composition prediction

The native composition owner implements the six fixed baselines from the
[held-out combination protocol](COMBINATIONS_PROTOCOL.md). It fits raw sparse
condition aggregates and predicts pairs of observed training targets in the same
declared context. It does not infer unseen-target effects, genetic interactions,
donor transfer, single-cell distributions or mechanistic parameters.

## Data and command boundaries

```text
source H5AD
  -> singlecell-h5ad-pseudobulk
  -> singlecell-composition-prepare
  -> training.json: control and declared single-target aggregates only
  -> singlecell-composition-fit
  -> singlecell-composition-predict: context and target identities only
```

Preparation reconstructs the native pseudobulk bundle with the current product
before selecting counts. Every selected condition must have exactly one aggregate,
and selected rows must agree on organism, cell group, donor and biological
replication identity. Single-target status is declared by the selection mapping;
the generic owner does not guess it from condition spelling. The Norman mapping
uses the 105 singles already audited in this release.

The training artifact records the original source/report fingerprints, evidence
classification, context, namespaces, target-to-condition mapping, feature IDs and
sparse counts. Matrix rows are **condition aggregates**, with control first and
then the declared singles. As in `VivoPseudobulkCounts`, the shared sparse storage's
`cellCount` field is its row dimension; no individual-cell identity is fabricated.

Fitting and prediction never receive paired-condition outcomes. Preparation may
read the full source to select training rows. Changing excluded conditions can
change source provenance fingerprints; it must not change selected counts or
fitted weights. Model verification reconstructs weights from the archived
training artifact. It does not retrieve experimental data by hash or authenticate
an independently supplied artifact's experimental claims.

```sh
numivivo singlecell-h5ad-pseudobulk original.h5ad --plan aggregation.json --output source-bundle
numivivo singlecell-composition-prepare source-bundle --plan selection.json --output training.json
numivivo singlecell-composition-fit training.json --output model
numivivo singlecell-composition-verify model
numivivo singlecell-composition-predict model --plan queries.json --output prediction
numivivo singlecell-composition-prediction-verify prediction
```

Destinations must not exist. Preparation writes a training file; fit and predict
write bundles. Prediction bundles contain the complete model reference, its
training file and receipt, the query plan, report and prediction receipt.
Verifiers use private snapshots, verify fingerprints and reexecute the relevant
calculation. Rehashing a changed model/report does not bypass reconstruction.

## Context and model contract

Selection schema version 1 contains `context`, `controlCondition`, `targets` and
`provenance`. Each target has an opaque `id` and source `condition`. Context has
`id`, `organism`, `featureNamespace`, `perturbationNamespace` and `countUnit`.
The query context must match every field exactly; its `queries` contain unique
query `id` values and two distinct observed target IDs. Unknown settings and
unknown/repeated targets are rejected.

Target IDs are separate from response-feature IDs. In Norman, target IDs are
the original guide-target names and features are the source Ensembl gene IDs.
This does not silently resolve the three unmatched target-name aliases identified
by the input audit. Measured single-condition responses can still be composed
without inventing that gene-identity mapping.

The frozen model retains control and single-condition CPM with complete source
feature denominators. Targets are sorted canonically before fitting. The owner
computes the fixed log/CPM composition rules and clips negative predicted
expression to zero. It preserves implied CPM totals and clipping diagnostics;
predictions are not reclosed or converted to counts.
Clipping counts describe floating-point operations, including tiny roundoff near
zero. In the mean-single reference, 226 features per query are clipped only
because of negative results no larger than 3.47e-17 in magnitude; this is not
evidence of biological suppression.

For each method, `expression[queryRows[i]]` gives the full gene vector for query i.
No-change and mean-single baselines share one stored row. Report feature/query
order is explicit. Unclipped quantities can be reconstructed from the frozen
control/target CPM and the versioned formulas; no held-out outcomes are needed.

## Resource limits

There are 2–256 training targets, 1–100,000 features, at most ten million dense
condition-feature values and five million sparse training nonzeros. Every
library must be positive, add without UInt64 overflow and total at most 2^53
for exact conversion to Double. Queries are limited to 1–256 pairs and twenty
million materialized prediction values across the shared and query-specific
methods. These bounds are checked before allocating those arrays.

Encoded training/model/query/report caps are 64/128/2/512 MiB. The CLI's common
plan reader additionally limits selection and query input to 128 KiB. Condition
CPM banks, predictions and JSON remain resident. This is not a million-cell or
GPU qualification; the input counts are obtained through the separately qualified
streaming H5AD route.

## Real-data qualification

`check_native_composition.py` accepts a fresh native Norman pseudobulk bundle and
the actual product binary. It verifies exact selected counts and feature order
against the published training archive, then checks model/prediction reconstruction,
all 131 pairs and 33,694 genes against the frozen external predictions, target
and training-order replay, invalid contexts/identities, overwrite and rehashed
model tampering. The external outcomes and baseline performance remain recorded
in [COMBINATIONS.md](COMBINATIONS.md).

On the full Norman release, all 786 query/method prediction arrays matched the
external reference exactly (maximum per-gene difference zero). Clipping counts
also matched exactly. CPM-total diagnostics differed by at most 6.8e-8 because of
summation order. All selected training counts and feature IDs matched exactly;
reordered training and reversed target order produced byte-identical models and
reports respectively. Native model and prediction reconstruction passed.

The release build and 54 single-cell tests in 15 suites passed. The full native
qualification exercised 20 commands, including 13 expected rejections. Additional
integrity checks exercise preparation overwrite and rehashed prediction tampering.
The main check includes rehashed model tampering; hashes alone are not treated as
proof of correct derived values.

```sh
python Tools/Omics/PerturbationPrediction/Norman/check_native_composition.py --binary /absolute/path/numivivo --pseudobulk source-bundle --out native-checks
python Tools/Omics/PerturbationPrediction/Norman/check_native_composition_integrity.py --binary /absolute/path/numivivo --pseudobulk source-bundle --qualified native-checks --out integrity-checks
```

Frozen native training/model/report artifacts and receipts are retained as gzip
files under `evidence/2026-09-09-native-composition`, alongside command logs,
numerical comparisons, build/test records and source hashes. Restore gzip files
to their original JSON names to reconstruct a bundle with the recorded product.
The model reference is stored once; the prediction's original reference copy was
byte-identical. The large source H5AD and source report are reused from the prior
full-source qualification by their exact hashes rather than duplicated here.
