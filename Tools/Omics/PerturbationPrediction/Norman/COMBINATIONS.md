# Norman held-out combination baseline evidence

All 131 paired-target conditions were withheld together. The predictors used
only control and the 105 single-target conditions, with all 33,694 source genes.
The paired outcomes were scored after training, prediction, frozen-model replay
and leakage checks completed. The [protocol](COMBINATIONS_PROTOCOL.md) fixes the
methods, transformations, panels and interpretation before scoring.

These are executable **external** baseline results and a reference for the next
native composition predictor. They do not implement that native predictor or
qualify unseen-target genes, genetic interactions, donors, cell types, Bayesian
uncertainty or mechanistic coupling.

## Results

Mean per-query response RMSE on natural-log(1 + CPM), with equal query weights:

| Baseline | All 33,694 genes | Training top 1,000 expression | Query top 200 training responses | Pairs worse than no change, all genes |
|---|---:|---:|---:|---:|
| No change | 0.192617 | 0.183185 | 0.974013 | 0 |
| Mean single response | 0.180105 | 0.152218 | 0.867338 | 15 |
| Additive log response | 0.187476 | **0.095324** | 0.604401 | 60 |
| Mean constituent log response | **0.162827** | 0.100689 | **0.517585** | 6 |
| Additive CPM response | 0.189615 | 0.110699 | 0.639785 | 59 |
| Shuffled-target additive log response | 0.240434 | 0.187618 | 0.836790 | 118 |

The training-expressed panel contains 19,429 genes and is also fully reported.
The query top-200 panel is selected from the additive responses of the two
already-observed singles; it uses no paired outcome. It highlights effects
expected under additivity, so it must be read alongside the other panels.

The additive-log baseline improves on its matched shuffled-target control in
121 of 131 pairs on the all-gene panel. The mean-constituent baseline has the
lowest mean all-gene error, but is worse than no change for
`BCL2L11_BAK1`, `BCL2L11_TGFBR2`, `KIF18B_KIF2C`, `PLK4_STIL`,
`PTPN12_OSR2` and `TMSB4X_BAK1`. Full log-additivity performs best on the
most-expressed panel but loses to no change for 60 pairs across all genes.
These failures are retained; no method was retuned after scoring.

Expression correlation alone would be misleading here: no change has mean
all-gene expression Pearson 0.991831 while its response correlation is undefined
because it predicts zero change. Additive log response has mean response
Pearson 0.525553; mean constituent response has 0.521679. Full per-query response
and expression metrics, including sign agreement and explained response sums
of squares, are in `scores/results.json`.

The predictions remain point estimates on a transformed expression scale.
For example, additive log predictions imply total CPM between 989,006 and
1,142,380; mean-constituent predictions imply 984,834 to 999,654. No reclosure or
conversion to counts is performed. Clipping counts and implied CPM sums are
retained for every method and pair in `predictions/diagnostics.json`.

The cell line and experimental pool are shared across conditions. Eight GEM
groups do not supply independent biological replication. This result therefore
does not provide biological-replicate confidence intervals or establish causal
interaction effects. The differing performance of fixed composition rules is
a measured baseline result, not a claim that realistic combinatorial biology
has been solved.

## Validation

* All 786 query-method prediction arrays reproduce exactly from the frozen model.
* Mutating every one of the 131 sealed paired outcomes leaves training arrays,
  model arrays and predictions unchanged.
* Reversing target order and training-target order preserves predictions.
* The independently fitted scikit-learn least-squares identity design agrees
  exactly with additive log predictions (maximum absolute difference zero).
* Ten negative cases reject unknown/repeated/wrong-count targets, negative,
  fractional, nonfinite and out-of-exact-range counts, and empty libraries.
* Independent complete prediction and scoring runs reproduce every saved array,
  receipt, diagnostic and metric exactly.

The raw count input is the already qualified, complete Norman condition
aggregation. Its reference archive SHA-256 is pinned by the splitter. Source
H5AD provenance, every-cell count/QC verification, and the source download
instructions remain in the [input qualification](README.md).

## Reproduce

Use the recorded NumPy/scikit-learn environment and fresh output directories.
The stage boundaries are deliberate: `predict` accepts no paired-outcome path.

```sh
python Tools/Omics/PerturbationPrediction/Norman/combinations.py split --reference Tools/Omics/PerturbationPrediction/Norman/evidence/2026-09-09/reference.npz --out inputs
python Tools/Omics/PerturbationPrediction/Norman/combinations.py predict --training inputs/training.npz --queries inputs/queries.json --out predictions
python Tools/Omics/PerturbationPrediction/Norman/check_combinations.py --reference Tools/Omics/PerturbationPrediction/Norman/evidence/2026-09-09/reference.npz --inputs inputs --predictions predictions --out checks
python Tools/Omics/PerturbationPrediction/Norman/combinations.py score --truth inputs/sealed-outcomes.npz --predictions predictions --out scores
```

Each prediction NPZ stores `expression` and `queryRows`; for query i, its full
per-gene prediction is `expression[queryRows[i]]`. Shared baselines store one
row rather than duplicating identical arrays. Feature order comes from the
frozen model. Unclipped quantities are reconstructed from that model's control
CPM and single-target responses through `predict_one`; paired outcomes are not
needed. Receipts bind inputs, implementation, protocol and output hashes.

Evidence is in `evidence/2026-09-09-combinations`. It includes the training-only
archive, sealed outcomes, frozen model, predictions, panels, diagnostics,
metrics, hashes and replay checks. No deep-learning model or published paper's
specific data processing was reproduced by this run.
