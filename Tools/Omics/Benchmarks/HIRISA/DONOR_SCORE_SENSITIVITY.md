# Sensitivity of completed prediction scores to individual donors

This retrospective diagnostic uses all **199 completed HIRISA folds**: 79 known-treatment folds and 120 cross/within-preparation folds. It reconstructs the published summaries from hash-verified score archives. No counts, models, folds or predictions change.

For each comparison, subtract baseline RMSE from candidate RMSE for each donor, then omit each donor score once and average the remaining differences. A negative difference favors the candidate. A stable improvement stays strictly negative in the full mean and every omission; stable worse stays strictly positive. Otherwise the result is donor-sensitive (or all-tied). These are score-aggregation sensitivities, not retrained leave-one-out experiments or confidence intervals.

| Ridge comparison setting | Baseline | Full-mean wins | Stable improvement | Stable worse | Donor-sensitive | Individual donor wins |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Known treatment | No change | 14/16 | 13 | 1 | 2 | 66/79 |
| Known treatment | Mean response | 4/16 | 3 | 10 | 3 | 31/79 |
| Cross preparation | No change | 12/12 | 12 | 0 | 0 | 56/60 |
| Cross preparation | Mean response | 3/12 | 3 | 8 | 1 | 21/60 |
| Matched within preparation | No change | 12/12 | 12 | 0 | 0 | 60/60 |
| Matched within preparation | Mean response | 7/12 | 6 | 4 | 2 | 38/60 |

For known treatment, the three stable ridge improvements over the mean are Monocyte IFNa, IFNb and IFNg. Its additional full-mean win, NK IFNa, changes sign under a donor omission. NK IFN-L1 is similarly sensitive against no-change. All twelve cross-preparation mean wins over no-change survive every omission, although four individual ridge predictions are worse. All twelve cross-versus-within ridge penalties also survive. This supports conditional average-response signal and a limited benefit from donor-context ridge, with the original failures retained.

All rows below use **18,082 source genes** and equally weighted donors. The JSON retains both original feature families, mean/median/ridge comparisons, every paired difference and every omitted-donor result.

| Setting / contrast | Donors | Ridge − no change | Omission range | Status | Ridge − mean | Omission range | Status |
| --- | ---: | ---: | --- | --- | ---: | --- | --- |
| Bcell / IFN-L1 | 5 | -0.078335 | [-0.080562, -0.073478] | stable-improvement | +0.001808 | [+0.001191, +0.002333] | stable-worse |
| Bcell / IFNa | 5 | -0.175803 | [-0.180085, -0.171393] | stable-improvement | +0.000639 | [+0.000463, +0.000787] | stable-worse |
| Bcell / IFNb | 5 | -0.191641 | [-0.195399, -0.187997] | stable-improvement | +0.000330 | [+0.000210, +0.000403] | stable-worse |
| Bcell / IFNg | 4 | -0.164541 | [-0.173570, -0.158630] | stable-improvement | +0.001207 | [+0.001054, +0.001567] | stable-worse |
| Monocyte / IFN-L1 | 5 | +0.011690 | [-0.002423, +0.022433] | donor-sensitive | +0.005107 | [+0.000782, +0.007692] | stable-worse |
| Monocyte / IFNa | 5 | -0.183504 | [-0.216833, -0.171173] | stable-improvement | -0.016787 | [-0.020943, -0.011219] | stable-improvement |
| Monocyte / IFNb | 5 | -0.196673 | [-0.241677, -0.182375] | stable-improvement | -0.014222 | [-0.017537, -0.009697] | stable-improvement |
| Monocyte / IFNg | 5 | -0.302409 | [-0.347825, -0.285423] | stable-improvement | -0.020802 | [-0.027912, -0.012007] | stable-improvement |
| NK / IFN-L1 | 5 | -0.000356 | [-0.002305, +0.002396] | donor-sensitive | +0.001642 | [+0.000629, +0.002575] | stable-worse |
| NK / IFNa | 5 | -0.174725 | [-0.182205, -0.169188] | stable-improvement | -0.000801 | [-0.001994, +0.000380] | donor-sensitive |
| NK / IFNb | 5 | -0.194380 | [-0.200894, -0.189740] | stable-improvement | +0.000129 | [-0.001148, +0.001101] | donor-sensitive |
| NK / IFNg | 5 | +0.008438 | [+0.005069, +0.010246] | stable-worse | +0.002984 | [+0.001642, +0.004112] | stable-worse |
| Tcell / IFN-L1 | 5 | -0.017720 | [-0.022385, -0.011981] | stable-improvement | +0.000211 | [-0.000424, +0.000620] | donor-sensitive |
| Tcell / IFNa | 5 | -0.131528 | [-0.136224, -0.126740] | stable-improvement | +0.003793 | [+0.001836, +0.004915] | stable-worse |
| Tcell / IFNb | 5 | -0.146243 | [-0.150541, -0.140651] | stable-improvement | +0.001326 | [+0.000916, +0.001852] | stable-worse |
| Tcell / IFNg | 5 | -0.025336 | [-0.034128, -0.021312] | stable-improvement | +0.000482 | [+0.000238, +0.000734] | stable-worse |
| cross / PBMC / B | 5 | -0.054308 | [-0.075238, -0.043654] | stable-improvement | -0.000933 | [-0.001101, -0.000661] | stable-improvement |
| cross / PBMC / CD4-T | 5 | -0.079852 | [-0.101907, -0.070759] | stable-improvement | +0.000147 | [-0.000242, +0.000738] | donor-sensitive |
| cross / PBMC / CD8-T | 5 | -0.122406 | [-0.134460, -0.117078] | stable-improvement | +0.006987 | [+0.006098, +0.008635] | stable-worse |
| cross / PBMC / Mono | 5 | -0.121384 | [-0.126948, -0.115082] | stable-improvement | +0.014317 | [+0.009367, +0.020347] | stable-worse |
| cross / PBMC / NK | 5 | -0.084507 | [-0.102348, -0.076829] | stable-improvement | +0.000494 | [+0.000266, +0.000747] | stable-worse |
| cross / PBMC / other-T | 5 | -0.044538 | [-0.055209, -0.039477] | stable-improvement | +0.001084 | [+0.000494, +0.001493] | stable-worse |
| cross / enriched / B | 5 | -0.084589 | [-0.085472, -0.083484] | stable-improvement | +0.002797 | [+0.001683, +0.003397] | stable-worse |
| cross / enriched / CD4-T | 5 | -0.069811 | [-0.076235, -0.065140] | stable-improvement | +0.003916 | [+0.001763, +0.005589] | stable-worse |
| cross / enriched / CD8-T | 5 | -0.030955 | [-0.039854, -0.020722] | stable-improvement | +0.004252 | [+0.003566, +0.005061] | stable-worse |
| cross / enriched / Mono | 5 | -0.099631 | [-0.121451, -0.082854] | stable-improvement | -0.014109 | [-0.018447, -0.007012] | stable-improvement |
| cross / enriched / NK | 5 | -0.085249 | [-0.094009, -0.078947] | stable-improvement | +0.005187 | [+0.004746, +0.005893] | stable-worse |
| cross / enriched / other-T | 5 | -0.041405 | [-0.044796, -0.039404] | stable-improvement | -0.002406 | [-0.002798, -0.002175] | stable-improvement |
| within / PBMC / B | 5 | -0.100846 | [-0.117887, -0.093526] | stable-improvement | -0.001025 | [-0.001569, -0.000003] | stable-improvement |
| within / PBMC / CD4-T | 5 | -0.126326 | [-0.155881, -0.115782] | stable-improvement | -0.004060 | [-0.005243, -0.002163] | stable-improvement |
| within / PBMC / CD8-T | 5 | -0.186306 | [-0.197012, -0.179551] | stable-improvement | -0.003165 | [-0.005111, +0.000092] | donor-sensitive |
| within / PBMC / Mono | 5 | -0.216530 | [-0.233112, -0.204253] | stable-improvement | -0.013167 | [-0.018506, -0.007712] | stable-improvement |
| within / PBMC / NK | 5 | -0.110468 | [-0.128674, -0.101868] | stable-improvement | -0.002515 | [-0.003245, -0.001365] | stable-improvement |
| within / PBMC / other-T | 5 | -0.064325 | [-0.075053, -0.058294] | stable-improvement | -0.007824 | [-0.008732, -0.006335] | stable-improvement |
| within / enriched / B | 5 | -0.183587 | [-0.188469, -0.178947] | stable-improvement | +0.001091 | [+0.000883, +0.001222] | stable-worse |
| within / enriched / CD4-T | 5 | -0.141076 | [-0.146503, -0.135198] | stable-improvement | +0.002162 | [+0.000361, +0.002994] | stable-worse |
| within / enriched / CD8-T | 5 | -0.118345 | [-0.125525, -0.112897] | stable-improvement | +0.004172 | [+0.001385, +0.005393] | stable-worse |
| within / enriched / Mono | 5 | -0.191191 | [-0.221927, -0.179194] | stable-improvement | -0.017686 | [-0.022008, -0.011481] | stable-improvement |
| within / enriched / NK | 5 | -0.182405 | [-0.191187, -0.177155] | stable-improvement | +0.000092 | [-0.001577, +0.001421] | donor-sensitive |
| within / enriched / other-T | 5 | -0.052023 | [-0.059501, -0.049321] | stable-improvement | +0.002582 | [+0.000694, +0.003856] | stable-worse |

## Transfer penalty

| Query preparation / lineage | Cross ridge − within ridge | Omission range | Status |
| --- | ---: | --- | --- |
| PBMC / B | +0.046538 | [+0.042650, +0.051694] | stable-worse |
| PBMC / CD4-T | +0.046474 | [+0.041956, +0.053974] | stable-worse |
| PBMC / CD8-T | +0.063900 | [+0.060610, +0.068380] | stable-worse |
| PBMC / Mono | +0.095146 | [+0.085423, +0.106164] | stable-worse |
| PBMC / NK | +0.025960 | [+0.023278, +0.029386] | stable-worse |
| PBMC / other-T | +0.019788 | [+0.017416, +0.022250] | stable-worse |
| enriched / B | +0.098998 | [+0.095463, +0.103462] | stable-worse |
| enriched / CD4-T | +0.071265 | [+0.065727, +0.078173] | stable-worse |
| enriched / CD8-T | +0.087390 | [+0.082961, +0.094467] | stable-worse |
| enriched / Mono | +0.091560 | [+0.084138, +0.101125] | stable-worse |
| enriched / NK | +0.097156 | [+0.094389, +0.100120] | stable-worse |
| enriched / other-T | +0.010618 | [+0.008388, +0.014705] | stable-worse |

## Interpretation and reproduction

A stable mean improvement can coexist with individual donor failures. Models share training donors; lineages and treatments reuse donors. Do not count 199 folds, genes or cells as independent biological replicates. Omitting a scoring donor does not remove that donor from other fitted models. This post hoc diagnostic cannot establish predictive interval coverage, repair preparation/batch confounding or validate a new study.

The [original known-treatment results](PREDICTION_RESULTS.md), [preparation-transfer results](CONTEXT_TRANSFER.md) and [biological assessment](../../../../Documentation/BiologicalPrediction.md) retain their scientific limits and negative findings.

From the repository root, using standard-library Python and new output paths:

```sh
python3 Tools/Omics/Benchmarks/HIRISA/prediction_sensitivity.py --out /tmp/hirisa-sensitivity.json --report /tmp/hirisa-sensitivity.md
```

The [evidence manifest](evidence/2026-09-11-prediction-sensitivity/manifest.json) binds the complete diagnostic and independent arithmetic check. Source archive identities are embedded in the diagnostic.
