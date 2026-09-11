# Donor-held-out learned reference benchmark

The subsequent [native shared-panel Kang–Ding transfer](CrossStudy/README.md)
checks every cell in both directions. Numerical agreement passes, but macro-F1
0.613/0.666 and zero megakaryocyte recall fail the declared coarse-family targets;
within-study reference results do not qualify general transfer.

This benchmark combines native H5AD projection and training-only native sparse
PCA with two external scikit-learn reference classifiers. It establishes a target
for [native frozen reference mapping](NATIVE.md). The original reference runner
independently reconstructs training centers from counts and verifies that frozen
projection reproduces native training scores. Native fitting now retains those
centers explicitly.

The [predeclared protocol](PROTOCOL.md) specifies the four Baron folds and fixed
parameters. All 8,569 cells and 20,125 source features are retained across the
folds; one donor is held out at a time. No query labels or counts enter HVG/PCA or
classifier fitting. Expression arrays stay sparse. Dense arrays are limited to
scores, loadings, classifier parameters and class probabilities.

Run with the benchmark Python environment and the qualified native product:

```sh
python run.py --source original.h5ad --plan stream-plan.json --binary /path/to/numivivo --out /new/run
python run.py --source original.h5ad --plan stream-plan.json --binary /path/to/numivivo --out /new/repeat
python check_replay.py /new/run /new/repeat --out replay.json
```

The plan must use the original Baron `X`, `native_sample`, `donor`, `cell_type`
and observation-index barcode axes. Use the preserved `source-plan.json` for the full donor mapping. Source IDs and metadata must not be silently remapped.

The subsequent [native balanced-logistic path](Logistic/README.md) now reproduces
all 8,569 external candidate labels with independent probability/objective checks.
Macro-F1 improves over native kNN in all four donors, while overall accuracy
falls in two and the original acinar/Schwann failures remain. The original
external benchmark below remains its comparison, not an independent validation.

## First real-data results

| Held-out donor | kNN accuracy | kNN macro F1 | Balanced logistic accuracy | Balanced logistic macro F1 |
|---|---:|---:|---:|---:|
| human1 | 0.962829 | 0.798132 | 0.976252 | 0.919378 |
| human2 | 0.986079 | 0.804712 | 0.979118 | 0.817137 |
| human3 | 0.964216 | 0.805484 | 0.976422 | 0.848056 |
| human4 | 0.953952 | 0.796431 | 0.947045 | 0.812805 |

Every source class occurs in each training fold. Uniform kNN misses every held-out
T cell in all four folds. Balanced logistic misses all three human2 acinar cells
and the single human3 Schwann cell. Thus overall accuracy does not establish
reliable rare-cell mapping. Neither model emitted fitting warnings in the first
runs. Class probabilities are uncalibrated, author labels are not experimental
proof of identity, and this single study cannot qualify cross-study transfer or
novel-label rejection. Human4 disease is confounded with donor.

The evidence retains per-label metrics/confusion matrices, predictions and
probabilities, frozen training means/loadings/scores/labels, logistic coefficients,
normalization target, source/binary hashes, plans, receipts and native reports.
The full original H5AD is shared with the earlier projection evidence and identified
by exact hash in the manifest. Projected files can be reconstructed from it.
Initial export used pandas object label arrays; the replay checker rejected them
because it disables pickle. The corrected runner explicitly saves Unicode labels.

API references: [Scanpy HVG](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html),
[scikit-learn logistic regression](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html),
[k-nearest neighbors](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html).

Corrected full-run replay passed all 32 comparisons: exact frozen/model arrays,
per-cell probabilities, metrics, native reports and projected H5AD bytes across
the four folds. Training HVG membership and PCA eigenvalue checks against
Scanpy passed in every fold. No Swift runtime source changed in the original benchmark commit `110e89d`.
