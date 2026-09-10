# HIRISA frozen donor-held-out prediction results

All 79 prespecified folds completed with native source reconstruction and byte-exact replay. Every model and query output was hashed before opening held-out treated counts for scoring. Independent NumPy reconstruction checks all selected features, control moments, coefficients and all four prediction vectors.

The maximum ridge treated-expression difference was 1.595e-13. Publication (including one full-source reconstruction and all fits) took 271.56 seconds; replay took 276.32 seconds on the physical M4 Pro Mac mini. Maximum RSS was 4,956,143,616/4,956,340,224 bytes. These are observed CPU timings, not a controlled scverse speed comparison or a Metal claim.

Mean, median and ridge responses improve the full-gene contrast RMSE over no-change in 14, 14 and 14 of 16 contrasts respectively. Ridge improves on the training-mean response in only 4 of 16. Monocyte IFN-L1 and NK IFNg retain worse errors than no-change for every learned baseline. No method is selected or promoted from these outcomes.

The tables report equally weighted held-out donor means within each contrast, in natural-log(1+CPM) units. Training has four donor pairs, except Bcell IFNg which has three. All 18,082 source genes and the separately frozen training-selected context family are shown. Both treated-expression and response RMSE are retained in the full results; subtracting the same control makes them equal apart from rounding. All fold-level MAE, response Pearson correlations and implied CPM sums remain in the archived JSON. Constant-response correlation, including no-change, remains null.

## All 18,082 source genes

| Population | Treatment | Donor folds | No change | Mean response | Median response | Context ridge |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Bcell | IFNa | 5 | 0.274165 | 0.097723 | 0.098303 | 0.098362 |
| Bcell | IFNb | 5 | 0.292331 | 0.100360 | 0.101751 | 0.100690 |
| Bcell | IFNg | 4 | 0.273697 | 0.107949 | 0.111332 | 0.109156 |
| Bcell | IFN-L1 | 5 | 0.178535 | 0.098393 | 0.098498 | 0.100201 |
| Monocyte | IFNa | 5 | 0.389989 | 0.223273 | 0.222991 | 0.206485 |
| Monocyte | IFNb | 5 | 0.421892 | 0.239441 | 0.237852 | 0.225219 |
| Monocyte | IFNg | 5 | 0.564236 | 0.282629 | 0.280472 | 0.261827 |
| Monocyte | IFN-L1 | 5 | 0.181475 | 0.188058 | 0.184872 | 0.193165 |
| NK | IFNa | 5 | 0.272105 | 0.098180 | 0.099425 | 0.097379 |
| NK | IFNb | 5 | 0.295204 | 0.100695 | 0.102043 | 0.100824 |
| NK | IFNg | 5 | 0.085704 | 0.091158 | 0.089881 | 0.094143 |
| NK | IFN-L1 | 5 | 0.087990 | 0.085992 | 0.087533 | 0.087634 |
| Tcell | IFNa | 5 | 0.235514 | 0.100192 | 0.099712 | 0.103986 |
| Tcell | IFNb | 5 | 0.246000 | 0.098431 | 0.099104 | 0.099757 |
| Tcell | IFNg | 5 | 0.131317 | 0.105499 | 0.105396 | 0.105981 |
| Tcell | IFN-L1 | 5 | 0.129138 | 0.111208 | 0.110893 | 0.111419 |

## Training-selected context genes

| Population | Treatment | Donor folds | No change | Mean response | Median response | Context ridge |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Bcell | IFNa | 5 | 0.292197 | 0.103741 | 0.104350 | 0.104421 |
| Bcell | IFNb | 5 | 0.312437 | 0.106823 | 0.108303 | 0.107171 |
| Bcell | IFNg | 4 | 0.299983 | 0.117342 | 0.121039 | 0.118664 |
| Bcell | IFN-L1 | 5 | 0.189982 | 0.104304 | 0.104413 | 0.106244 |
| Monocyte | IFNa | 5 | 0.436558 | 0.249736 | 0.249424 | 0.230947 |
| Monocyte | IFNb | 5 | 0.472964 | 0.268247 | 0.266468 | 0.252299 |
| Monocyte | IFNg | 5 | 0.632900 | 0.316814 | 0.314379 | 0.293476 |
| Monocyte | IFN-L1 | 5 | 0.202680 | 0.210065 | 0.206481 | 0.215771 |
| NK | IFNa | 5 | 0.308489 | 0.110694 | 0.112115 | 0.109781 |
| NK | IFNb | 5 | 0.335324 | 0.113680 | 0.115215 | 0.113828 |
| NK | IFNg | 5 | 0.096503 | 0.102696 | 0.101255 | 0.106085 |
| NK | IFN-L1 | 5 | 0.099203 | 0.096905 | 0.098678 | 0.098757 |
| Tcell | IFNa | 5 | 0.264906 | 0.112082 | 0.111536 | 0.116381 |
| Tcell | IFNb | 5 | 0.276285 | 0.109252 | 0.110008 | 0.110753 |
| Tcell | IFNg | 5 | 0.146237 | 0.117350 | 0.117244 | 0.117888 |
| Tcell | IFN-L1 | 5 | 0.143682 | 0.123614 | 0.123279 | 0.123842 |

## Evidence and limits

The [prediction archive](evidence/2026-09-10-prediction) retains every native model/prediction, source and implementation bindings, original protocols/folds, both output freezes, numerical checks, all scores, product tests and failed attempts. The first transport preparation missed the pool component of the batch identity and was repaired before publication. A later transport description contained the wrong GEO accession; that native attempt was explicitly interrupted during source reconstruction before fold fitting. Its records remain alongside the corrected GSE306664 run. Neither repair changed the frozen scientific protocol, folds or numerical model.

This is empirical prediction of a known perturbation in an unseen donor from its observed control, using experimental enrichment populations rather than verified single-cell annotations. The separate PBMC transfer split is still unfrozen. The results do not establish unseen perturbation identity, calibrated uncertainty, causal mechanisms or general biological validity. Metadata/report residency and million-cell PCA/graph/integration remain open.

Score JSON SHA256: `8d49cd4eb3ce2ced43f006cec931d848903251cfff0d496621004cc3637ce44f`.

Reproduce the numerical checks and scores with `python score_prediction_batch.py --root /absolute/path/hirisa`, then render this report with `python report_prediction.py --root /absolute/path/hirisa --output PREDICTION_RESULTS.md`. Restore the complete native source/report and frozen bundle according to the archive links; verification of a historical native receipt requires its original executable identity.
