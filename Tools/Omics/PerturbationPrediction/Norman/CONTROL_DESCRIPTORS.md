# Control-only gene descriptors: second transfer failure

Untreated-cell correlation descriptors still fail to improve on the mean-single
response baseline for unseen targets in the full Norman benchmark. Their all-gene
average is also slightly worse than the matched shuffled descriptor control.
This experiment does not support promoting this model as a useful native
unseen-target predictor.

The [protocol](CONTROL_DESCRIPTORS_PROTOCOL.md) was specified after the first
co-response failure and before scoring this model. These results reuse the
research benchmark, so they are not an untouched external validation. No
hyperparameter or descriptor selection was changed after scoring.

## Inputs and coverage

The source is the same byte-pinned Norman H5AD and independently native-qualified
condition reference as the [full input qualification](README.md). Descriptors
use only 11,855 control cells. All 33,694 genes are streamed for control moments;
2,000 landmarks are chosen by control log1p-CPM variance, requiring at least ten
expressing controls and excluding all directly matched intervention target genes.
This is a declared variance selector, not a Scanpy HVG flavor. Sparse selected
control columns contain 12,988,400 stored entries. Cross-products produce a small
target-by-landmark correlation matrix; no dense cells-by-genes matrix is built.

97 targets have unique exact symbol matches and sufficient control expression.
The three unresolved aliases remain C19orf26, C3orf72 and KIAA1804. Five other
targets fail the predeclared ten-control-cell threshold:

| Target | Expressing control cells |
| --- | ---: |
| FOXF1 | 8 |
| HNF4A | 7 |
| HOXA13 | 3 |
| MAP2K6 | 5 |
| PTPN13 | 0 |

Those eight targets are unsupported for descriptor predictions, not replaced
with an unknown-target fallback. Both generic baselines still score all 105.
Each fold excludes every single/pair involving its target; fitting uses only
control and other singles. Control-cell readouts of the held gene are available,
but its perturbation response is never passed to the fitter.

## Results

Mean per-target response RMSE, matched on the same 97 supported targets:

| Method | All genes | Training top 1,000 | Worse than no-change, all genes |
| --- | ---: | ---: | ---: |
| No change | 0.132701 | 0.122199 | 0 |
| Mean training-single response | 0.126328 | 0.105049 | 22 |
| Control-correlation ridge | 0.129260 | 0.107868 | 29 |
| Shuffled control-correlation ridge | 0.129039 | 0.108210 | 28 |

Ridge beats the generic mean baseline on only 19/97 all-gene targets and 37/97
training-top-1,000 targets. It beats the matched shuffled model on 44/97 and
47/97, respectively. The shuffle is a fixed rotation, not a significance test.
The generic mean uses all 104 available singles; ridge and its matched shuffled
control use the 96 other descriptor-supported training targets. Their training
response means therefore need not be identical. Scores, correlations, sign
agreement, clipping counts and implied CPM sums retain every individual failure.

The earlier descriptor experiment had 102 supported targets. Comparing its
published averages directly with these 97-target averages would mix cohorts.
Removing its self-intervention mismatch has not demonstrated better predictive
transfer. Control correlations remain observational and can encode cell-state
or technical confounding; these results do not establish regulatory causality.

## Verification

Two complete preparations produced byte-identical descriptors, coverage and
receipts. Every control cell's reconstructed total matches the pinned reference;
barcodes, feature identities and condition labels match the original H5AD.
Three independent two-vector NumPy correlation checks have maximum absolute
error 7.78e-16. Full prediction runs reproduce arrays, fold records and receipts
byte for byte. Dual ridge and independent scikit-learn SVD Ridge agree within
1.45e-14 across all supported folds, for ordinary and shuffled outputs.

Changing every excluded outcome leaves all 105 fitting inputs identical. Full
refits on AHR, KIAA1804 and ZNF318 reproduce their predictions and
panels exactly after those mutations; KIAA1804 exercises the unsupported-target
generic baselines, while AHR and ZNF318 also exercise descriptor ridge. No landmark is an intervention target.
Unsupported rows are explicit NaNs and are excluded from descriptor scoring.
Both complete scoring runs are identical.

The predictor's orchestrator reads the pinned reference to select folds; only
selected training data enters fitting. This is function-level outcome isolation,
not a separate operating-system sandbox. Descriptors and prediction outputs are
frozen before the scoring invocation. Prediction artifacts are not reusable
native model bundles. Swift code and native product qualification are unchanged.

Observed local Apple M4 timings: descriptor preparation 26.48 seconds, peak RSS
1,010,122,752 bytes; first complete prediction 27.67 seconds, peak RSS
2,011,906,048 bytes. Independent SVD checks are included. These are observed run
costs; some independent preparation/prediction work overlapped, and these are
not isolated speed comparisons. Packages and all logs are retained.

## Reproduce

From this directory in the established benchmark Python environment:

```sh
python control_descriptors.py prepare --source /path/original.h5ad --reference /path/reference.npz --out descriptors
python control_descriptors.py predict --reference /path/reference.npz --descriptors descriptors --out predictions
python control_descriptors.py score --reference /path/reference.npz --predictions predictions --out scores
```

Repeat all three commands with new output destinations, then run:

```sh
python check_control_descriptors.py --reference /path/reference.npz --descriptors descriptors --repeat-descriptors repeat-descriptors --first predictions --repeat repeat-predictions --out checks.json
```

Evidence in `evidence/2026-09-09-control-descriptors` retains descriptors,
coverage, predictions, all per-target scores, checks, receipts and logs. Original
H5AD/reference files remain pinned by existing qualifications, without duplication.
Better target information, independent context validation, native unseen-target
support and Bayesian/mechanistic coupling remain open.
