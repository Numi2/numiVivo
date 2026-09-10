# Program decoder calibration: post-result development

**PARTIAL: 28/32 controls are sensitive and all four integration candidates meet
the original loss margins under a within-library fit. Four controls remain
insufficient, so the complete preservation gate still does not pass.** This is
a separately declared development experiment after inspecting the [original
mixed-objective failures](INTEGRATION_PROGRAMS.md). Those results remain intact;
this follow-up is not independent biological replication or a replacement gate.

## What changed and why

The original ridge decoder fits total RNA program variation, including large
between-library means. Its primary score then centers both truth and prediction
within each held-out library. The training and evaluation objectives therefore
emphasize different variation. Negative original baseline R² and sensitivity to
integration can reflect that mismatch, not only lost program information.

This follow-up centers training PCs and measured RNA separately within each
training library, then fits equal-library-weighted ridge on that variation.
PC scaling uses pooled within-library training variance. Ridge remains 1; all
20 PCs, programs, feature coverage, original donor splits, margins and candidate
matrices remain unchanged. The global intercept uses training means only and
serves secondary uncentered diagnostics. Held-out target means enter scoring
only. Full-cohort PCA/integration still make this a transductive information
preservation diagnostic, not prospective treated-outcome prediction.

Protocol SHA256 is
`e27f02aec5ed1d00ce03eeddcba6ddecb44dca5f3c91bac9928f2c1f57de34c0`.
It was frozen before this follow-up's metrics, but after the original results
were known. No new score chose a different regularizer, threshold, program,
cohort or winner. Every original cell remains represented in the stored moments:
1,612,594 total cells, 131 libraries; 79 matched folds score 1,275,710 distinct
cells, with the other 336,884 accounted for outside these folds.

The evaluator reuses exact previously verified full-source sufficient statistics
and first reproduces every original mixed-objective coefficient/intercept exactly.
It does not rerun RNA preparation or integration. New full-cell checks independently
read the original latent matrices and frozen RNA targets, directly center actual
training rows, and solve augmented weighted least squares using SVD. Thus the
numerical qualification does not rely on the reused moments alone.

## Results and remaining insufficient controls

The original objective had 14/32 sensitive controls and three sensitive margin
failures for every candidate. The aligned objective has 28/32 sensitive controls
and no mean or individual-fold margin failure in any candidate. Baseline mean
within-library R² now ranges from **0.025845 to 0.502394**; the four insufficient
comparisons are the five-gene program in NK IFNg, NK IFN-L1, Tcell IFNg and Tcell
IFN-L1, with baseline skill 0.027189, 0.025845, 0.034887 and 0.041149 respectively.
Their failure to exceed the unchanged 0.05 control threshold remains explicit.

| Candidate | Maximum contrast-mean R² loss | Maximum individual-fold R² loss | Complete gate |
| --- | ---: | ---: | --- |
| native | 0.006099 | 0.031074 | Insufficient controls |
| harmony-7 | 0.006683 | 0.035727 | Insufficient controls |
| harmony-19 | 0.006945 | 0.035695 | Insufficient controls |
| harmony-41 | 0.006084 | 0.032640 | Insufficient controls |

Allowed losses remain 0.05 mean and 0.10 per fold. Identity is exactly unchanged.
Library-mean-only erasure removes all within-library signal and is detected in
all 28 sensitive comparisons. Nine control tests pass: three new centered-fit,
training-offset-invariance and held-out-independence tests plus the six original
program controls. All five complete actual-data SVD checks pass; maximum R²
error is 3.55e-13 and maximum coefficient error is 6.40e-14. Every original fold
and both programs are checked, including direct held-out residuals.

The reused-statistic development evaluation took 0.332 seconds. Independent
per-cell checks each took 12.6–12.8 seconds with at most 798,064,640 resident
bytes. These are shared-host observations; cached fitting time excludes the
original RNA/PCA/integration/moment preparation and is not a native-method speed
comparison. The native candidate is the original seed-7 integration executable
`e360645362ac0707cb497eaa64277e5bf792b28d79bddea707421fdcd043495d`;
Harmony seeds 7/19/41 are the existing full runs. No new integration fit is claimed.

## All comparisons

H is the 93 matched Hallmark members; F is the overlapping five-gene program,
not independent replication. Each candidate column is mean / maximum-fold R²
loss. Original and calibrated baseline columns expose the fitting-objective
change without hiding weak or negative original results.

| Population | Treatment | Program | Original baseline R² | Calibrated baseline R² | Sensitive | Native 7 | Harmony 7 | Harmony 19 | Harmony 41 |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| Bcell | IFNa | H | -0.282533 | 0.218581 | yes | 0.002135 / 0.004769 | 0.001898 / 0.004179 | 0.001549 / 0.004029 | 0.001914 / 0.004470 |
| Bcell | IFNa | F | -0.609161 | 0.094383 | yes | 0.002635 / 0.006232 | 0.002593 / 0.005745 | 0.001030 / 0.002514 | 0.002782 / 0.006418 |
| Bcell | IFNb | H | -0.292730 | 0.228872 | yes | 0.002079 / 0.005091 | 0.001899 / 0.004638 | 0.001827 / 0.004732 | 0.001864 / 0.004852 |
| Bcell | IFNb | F | -0.613443 | 0.112439 | yes | 0.002555 / 0.005987 | 0.002551 / 0.005522 | 0.001186 / 0.003148 | 0.002744 / 0.006259 |
| Bcell | IFNg | H | 0.063957 | 0.217543 | yes | 0.002597 / 0.005412 | 0.002196 / 0.004625 | 0.002972 / 0.005450 | 0.002079 / 0.005005 |
| Bcell | IFNg | F | 0.045481 | 0.060031 | yes | 0.002157 / 0.005113 | 0.002001 / 0.004530 | 0.000673 / 0.001693 | 0.002137 / 0.004954 |
| Bcell | IFN-L1 | H | 0.132824 | 0.192141 | yes | 0.003369 / 0.009614 | 0.003684 / 0.009476 | 0.002774 / 0.009019 | 0.003526 / 0.009377 |
| Bcell | IFN-L1 | F | -0.005252 | 0.071064 | yes | 0.002788 / 0.006250 | 0.003489 / 0.006827 | 0.000926 / 0.006785 | 0.003539 / 0.006464 |
| Monocyte | IFNa | H | 0.461946 | 0.462832 | yes | -0.000098 / 0.025472 | -0.000101 / 0.028013 | -0.000021 / 0.029009 | -0.000017 / 0.026740 |
| Monocyte | IFNa | F | 0.355779 | 0.404661 | yes | 0.000350 / 0.014457 | 0.000826 / 0.017730 | 0.000835 / 0.018046 | 0.000762 / 0.015509 |
| Monocyte | IFNb | H | 0.463483 | 0.466565 | yes | -0.000016 / 0.024298 | 0.000122 / 0.027243 | 0.000146 / 0.028175 | 0.000167 / 0.026015 |
| Monocyte | IFNb | F | 0.345129 | 0.402604 | yes | 0.001289 / 0.016215 | 0.001700 / 0.019541 | 0.001775 / 0.019644 | 0.001692 / 0.017333 |
| Monocyte | IFNg | H | 0.415722 | 0.443062 | yes | -0.003552 / 0.019310 | -0.003181 / 0.023005 | -0.002721 / 0.024633 | -0.002940 / 0.021234 |
| Monocyte | IFNg | F | 0.326662 | 0.336592 | yes | -0.004077 / 0.007099 | -0.003161 / 0.012286 | -0.003136 / 0.010852 | -0.003141 / 0.008552 |
| Monocyte | IFN-L1 | H | 0.537780 | 0.502394 | yes | 0.001718 / 0.031074 | 0.001750 / 0.035727 | 0.002345 / 0.035695 | 0.001884 / 0.032640 |
| Monocyte | IFN-L1 | F | 0.455269 | 0.440016 | yes | 0.006099 / 0.022940 | 0.006683 / 0.028141 | 0.006945 / 0.026906 | 0.006084 / 0.023362 |
| NK | IFNa | H | -0.441691 | 0.147982 | yes | -0.001109 / 0.000586 | -0.000803 / 0.001194 | -0.001139 / 0.001054 | -0.001112 / 0.001090 |
| NK | IFNa | F | -0.959054 | 0.079258 | yes | -0.000598 / 0.001259 | -0.000484 / 0.001306 | -0.000656 / 0.001102 | -0.000632 / 0.001387 |
| NK | IFNb | H | -0.448468 | 0.165887 | yes | -0.000770 / 0.001413 | -0.000478 / 0.002077 | -0.000785 / 0.001879 | -0.000783 / 0.001899 |
| NK | IFNb | F | -0.983298 | 0.100245 | yes | -0.000451 / 0.001410 | -0.000345 / 0.001741 | -0.000533 / 0.001420 | -0.000507 / 0.001483 |
| NK | IFNg | H | 0.129715 | 0.115415 | yes | 0.000744 / 0.001852 | 0.001070 / 0.001836 | 0.000686 / 0.001277 | 0.000733 / 0.001809 |
| NK | IFNg | F | 0.028259 | 0.027189 | no | -0.000265 / 0.000390 | -0.000201 / 0.000397 | -0.000300 / 0.000321 | -0.000234 / 0.000432 |
| NK | IFN-L1 | H | 0.104062 | 0.094510 | yes | 0.001509 / 0.008855 | 0.001815 / 0.009069 | 0.001440 / 0.008351 | 0.001688 / 0.009325 |
| NK | IFN-L1 | F | 0.028581 | 0.025845 | no | 0.001277 / 0.006927 | 0.001302 / 0.006771 | 0.001084 / 0.006418 | 0.001317 / 0.007069 |
| Tcell | IFNa | H | -0.181607 | 0.142276 | yes | -0.000080 / 0.001686 | -0.000029 / 0.001902 | 0.000130 / 0.001765 | -0.000102 / 0.001962 |
| Tcell | IFNa | F | -0.560828 | 0.082251 | yes | -0.000101 / 0.000929 | -0.000071 / 0.001120 | -0.000075 / 0.000979 | -0.000121 / 0.001065 |
| Tcell | IFNb | H | -0.216580 | 0.152996 | yes | 0.000405 / 0.002292 | 0.000497 / 0.002620 | 0.000606 / 0.002321 | 0.000394 / 0.002670 |
| Tcell | IFNb | F | -0.638252 | 0.081538 | yes | -0.000014 / 0.000846 | 0.000025 / 0.000942 | 0.000010 / 0.001010 | -0.000037 / 0.000842 |
| Tcell | IFNg | H | 0.283946 | 0.269836 | yes | 0.001494 / 0.003142 | 0.001664 / 0.002986 | 0.001663 / 0.003776 | 0.001691 / 0.003336 |
| Tcell | IFNg | F | 0.037452 | 0.034887 | no | 0.000045 / 0.000682 | 0.000110 / 0.000585 | 0.000096 / 0.000736 | 0.000043 / 0.000681 |
| Tcell | IFN-L1 | H | 0.163441 | 0.158741 | yes | 0.001155 / 0.004832 | 0.001567 / 0.004796 | 0.001076 / 0.004673 | 0.001445 / 0.004895 |
| Tcell | IFN-L1 | F | 0.042632 | 0.041149 | no | 0.000271 / 0.005315 | 0.000315 / 0.005522 | 0.000157 / 0.004664 | 0.000287 / 0.005372 |

## Interpretation and evidence

These results show that the earlier losses are decoder-objective dependent.
They support retained linear program information in the 28 control-sensitive
comparisons under this development fit. They do not prove that all biological
information is preserved, that integration removes only technical variation,
or that the earlier decoder diagnosed causal biological damage. Four insufficient
controls, rare-cell/neighborhood preservation, native multi-seed robustness and
independent-context replication remain open. A new claim requires an appropriate
measured endpoint; neither a larger donor-classification score nor a renamed
gate can supply that evidence.

The separate count-based donor-response prediction results and [native RNA
program scoring](NATIVE_PROGRAM_RESULTS.md) are unchanged. Supplied program
scores remain descriptive expression summaries, not calibrated activities.

The [archive manifest](evidence/2026-09-10-program-calibration/manifest.json)
retains 37 members: frozen protocol and source identities, exact scripts/tests,
all seven results and complete fold/moment records, five SVD checks and execution
logs. It maps each reused input to its exact original archived result, including
the original failures. Large source/latent/RNA arrays remain in the earlier
study storage and archives under unchanged hashes.

```sh
python Tools/Omics/Benchmarks/HIRISA/verify_archive.py \
  Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-10-program-calibration
```
