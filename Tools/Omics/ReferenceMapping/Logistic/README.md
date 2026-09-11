# Native balanced logistic reference mapping

The existing reference workflow now optionally fits class-balanced multinomial
logistic regression in frozen training PCA coordinates. It learns the supplied
reference labels and emits candidate labels plus normalized, **uncalibrated**
class probabilities. Uniform kNN remains the default.

Add `logistic` to an existing `singlecell-reference-fit` plan:

```json
"logistic": {
  "penalty": 1,
  "gradientTolerance": 1e-8,
  "maximumIterations": 2000,
  "maximumWork": 20000000000
}
```

Use the existing commands to fit, map and reconstruct:

```sh
numivivo singlecell-reference-fit training.h5ad --plan fit.json --output reference
numivivo singlecell-reference-verify reference
numivivo singlecell-reference-map query.h5ad --plan query.json \
  --reference reference --output mapped
numivivo singlecell-reference-map-verify mapped
```

The [reference input contract](../NATIVE.md) still applies: explicit label/source
provenance, matching complete gene universe and organism, training-only sparse
HVG/PCA, query library normalization, and no query-label mapping. No novel label,
missing gene or unseen tissue is inferred. Empty query libraries remain unmapped.

## Model and numerical contract

Each training PCA component is centered and scaled using the unweighted training
mean and population variance. Numerically constant components use scale one.
For n training cells, K classes and n_c cells in class c, each cell receives
weight n/(K n_c). The convex objective is weighted mean negative log likelihood
plus `penalty/(2n)` times the squared coefficient norm. Intercepts are unpenalized.
The multinomial formulation uses one parameter vector per class, including for
two-class inputs. For the multi-class Baron comparison, penalty 1 matches the
original scikit-learn C=1 model.

A native L-BFGS solver uses ten history pairs and an Armijo line search. Only a
fit meeting the declared maximum absolute gradient is published. Work/iteration
limits and line-search failure are errors. The model retains training scale,
class counts/weights, parameters, accepted objective history, final gradient,
iterations, evaluations and charged work. The default gradient tolerance is
1e-7; the qualified example above requests 1e-8. Default penalty/iteration/work
limits are 1/2000/20 billion. Query plans may optionally set
`maximumClassifierOperations` (default two billion, maximum twenty billion).

Training work charges two `n*K*(components+1)` terms per objective/gradient
execution; query classification charges `K*(components+1)` per nonempty row.
These are explicit arithmetic budgets, not wall-time or FLOP measurements.
Existing PCA, projection and output limits remain independent. Dense arrays are
limited to cells-by-components, parameters/history and class outputs, with at
most 100,000 reference cells, 64 components and 256 classes. This is the resident
reference path, not million-cell or GPU qualification.

Logistic query rows contain `classProbabilities` in the report's class order;
neighbor indices, distances and vote arrays are empty because no neighbor search
was performed. Ties choose the first sorted class. `distanceOperations` is zero,
and `classifierOperations` is explicit. Omitting `logistic` preserves historical
model/report encoding: these optional logistic fields remain absent. The native
verifier refits from retained original counts/labels and reconstructs predictions;
changing model coefficients and rehashing the receipt does not bypass it.

Class balancing changes the fitted class prior. The probabilities are not
calibrated confidence for the query population and do not constitute identity
proof. Unknown-class rejection, independent-study validation and reliable rare
class mapping remain separate requirements.

## Complete donor-held-out qualification

The [protocol](PROTOCOL.md) retains all four original Baron donor splits, all
8,569 query cells and all 20,125 source features. Previous external classifier
results were already inspected; this is native implementation and requalification
on development data, not a new independent biological experiment. Human4 disease
status is confounded with donor. Source author labels are the evaluation target,
not experimentally re-established cell identities.

The first complete native run used gradient tolerance 1e-7 and passed all native
reconstructions. Its first independent probability check failed the fixed 1e-4
agreement bound (maximum difference 0.00014895). Outcome-label scoring had not
started. The [numerical refinement](NUMERICAL_REFINEMENT.md) preserves that run
and refits all four folds at 1e-8 with the same model, penalty and source data;
it does not loosen the comparison tolerance or tune to biological accuracy.

The refined run passes all **17 native fit/map/reconstruction/repeat commands**.
Six Swift tests and fifteen lifecycle/regression commands pass. All native
probabilities agree with explicit independent softmax within 1e-12 and with the
tighter scikit-learn fits within **2.24605e-5**, below the unchanged 1e-4 bound.
Every predicted label matches both the tighter fit and the earlier external
balanced-logistic benchmark. Training scale, class weights, objective and gradient
checks pass; repeated scoring reproduces all nine output files byte-for-byte.

| Held-out donor | Cells | kNN accuracy | Native logistic accuracy | kNN macro-F1 | Native logistic macro-F1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| human1 | 1,937 | 0.962829 | 0.976252 | 0.798132 | 0.919378 |
| human2 | 1,724 | 0.986079 | 0.979118 | 0.804712 | 0.817137 |
| human3 | 3,605 | 0.964216 | 0.976422 | 0.805484 | 0.848056 |
| human4 | 1,303 | 0.953952 | 0.947045 | 0.796431 | 0.812805 |

Macro-F1 improves in every donor, while overall accuracy decreases in human2
and human4. The model identifies six of seven held-out T cells, compared with
zero for kNN, but that tiny stratum is insufficient rare-cell validation. It still
misses all three human2 acinar cells and the single human3 Schwann cell. Complete
per-label precision/recall and confusion matrices are retained, including false
positives. No default classifier is replaced or general annotation claim promoted.

The actual scoped release executable is SHA-256
`2c0b07146cdbfe3c4d0e525a3354bca232a0214b5dfaa7875d363f5dbe8b01a0`,
using the previously qualified HDF5 library
`a00ffbf8ab94ad81f67231a1ae01df748689e1c35a3615a57f88f4710b8d213e`.
The refined input freeze is
`5ea4ea1c0005d23112fc692e9fdf43a5f0e7e1a9215e4f08e7625d4c0c1fd70a`
and prediction freeze is
`fd31a9835c5e1733c7f84bba4a4da410c70e91e8a65ec4a00227c56f6938fd76`.
Execution used the physical M4 Pro Mac mini with 24 GB, macOS 26.6 and Swift 6.3.3.
Maximum resident memory across the refined product commands was 248,545,280
bytes; this is an operational observation, not a comparative speed claim or
full-package qualification.

Retained failures include constant-column scaling detected by a Swift test,
the initial probability-accuracy mismatch, and a lifecycle checker that wrongly
expected a one-row empty fixture. The actual retained fixture contains one empty
and two nonempty cells; the repaired checker verifies both behaviors and exact
classification-work accounting. Native reconstruction rejects modified model
coefficients even after a receipt is rehashed. Full historical default kNN
model/report bytes remain identical.


## Reproduce and inspect evidence

Build the actual scoped product CLI and test the owner:

```sh
export NUMIVIVO_HDF5_LIBRARY=/absolute/path/libhdf5.dylib
bash Tools/Omics/H5AD/build.sh /path/to/runtime --with-cli
bash Tools/Omics/ReferenceMapping/Logistic/test.sh /path/to/runtime
```

`prepare.py` verifies and compresses the existing complete source-linked donor
splits, then freezes all original plans. `refine.py` creates a separate input
freeze for the stricter convergence setting, referencing the same gzip sources.
`run.py` executes every fit/map and native verifier before freezing predictions.
`score.py` then compares every probability with explicit softmax and a tighter
independent scikit-learn fit, checks the objective/gradient, and scores all author
labels. Its arguments select the original source-linked benchmark directory;
no dataset is downloaded or reduced by these recipes. `regression.py` exercises
restored full native bundles and the retained historical structural queries.

Full bundles are content-deduplicated gzip tar archives. `bundles.py ARCHIVE`
verifies all member identities; `bundles.py ARCHIVE --restore NEW_DIRECTORY`
reconstructs every file. On APFS, identical files use separate copy-on-write
clones, not shared mutable hardlinks. Native verification requires the recorded
executable. Both the first and refined numerical runs are retained. Finished
scratch is removed only after every archive object and reconstruction path is
verified and no process holds it open. Original sources remain intact.

[Compact evidence](evidence/2026-09-11/manifest.json) contains plans, source and
executable identities, complete numerical/label metrics, failures, native logs
and executed recipes. The full original and refined archives remain on both
hosts under `numivivo-reference-logistic-20260911/{native,native-qualified}`.
They preserve 5,175,394,884 logical file bytes in 333,868,228 compressed bytes,
including duplicate reference/source paths. APFS clones already shared some
physical storage, so this is a representation-size comparison, not a claim of
that much reclaimed disk. Completed build intermediates totaling 100,196,510
bytes were separately removed after structural and open-handle checks; all
qualified executables, source manifests and research inputs remain intact.
