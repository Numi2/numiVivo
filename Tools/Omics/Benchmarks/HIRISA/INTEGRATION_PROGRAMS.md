# Complete HIRISA program-preservation diagnostic

**FAIL on the frozen complete preservation gate for native seed 7 and all three
Harmony references.** All four candidates lose more than the declared margins
in the same three control-sensitive comparisons. Eighteen of 32 comparisons
have insufficient erasure-control sensitivity. The previous [coarse response
checks](INTEGRATION_RESPONSE.md) still pass; they did not test these gradients.
The separate [raw-count donor-response predictions](PREDICTION_RESULTS.md) are
unchanged and do not use integrated PCA.

## Cohort, timing and measured targets

All **1,612,594 cells, 18,082 genes and 131 libraries** remain in preparation and
library moments. The original 79 matched donor folds evaluate 1,275,710 distinct
cells; the other 336,884 are accounted for outside these folds. No gene-count
matrix is densified. Original source SHA256 is
`0873e698ebf8770a54e6dba09724ffbeda5e1a67dbf24d4e223552d4aca7969c`.

The two existing Kang program definitions are unchanged. Exact unique original
symbol-to-Ensembl mapping retains 93/97 Hallmark interferon-alpha members and
5/5 five-gene members. CASP1, HLA-C, PNPT1 and WARS1 are missing from the first
mapping and retained as missing. Coverage thresholds are 0.8 and 1.0. There is
no alias expansion or outcome-selected membership. The five-gene set overlaps
the Hallmark set and is not independent replication. These weighted expression
summaries are not calibrated pathway activities or authoritative cell labels.

Protocol SHA256 is
`90e55cab09b15de892208a21b43ca6cb8ddfd6811c5c5fac723645b4e518c1d9`.
It was frozen after prior coarse HIRISA outcomes, Kang failures and metadata-only
symbol coverage were known, but before HIRISA per-cell program scores or gradient
metrics were inspected. The archive retains this disclosure and execution
freezes. No new result selected a different margin, cohort or model.

Program scores use signed L1-weighted log1p expression normalized to 10,000,
with every source count contributing to each cell's total. The independent
stream reads all 3,845,991,249 stored entries and 7,700,096,227 UMIs, compares
every cell total with the prior independent source QC, and compares each score
with a separate SciPy sparse matrix product. Maximum discrepancy is 2.67e-15.
Preparation took 56.65 seconds and 218,169,344 bytes peak RSS. All cells have
nonzero totals. Missing/empty targets would remain unavailable rather than zero.

## Frozen diagnostic

For each representation and program, equal-library-weighted linear ridge uses
all 20 PCs, training-only centering/scaling, ridge 1 and an unpenalized intercept.
Each donor fold scores original RNA targets withheld from fitting. Full-cohort
PCA and integration precede the diagnostic, so it is transductive information
preservation, not prospective donor-response prediction.

Within-library R² is one minus centered prediction-error variance divided by
measured target variance. Both prediction and target are centered within each
held-out library at scoring. Fitting uses the full weighted training response,
including between-library differences. Frozen allowed losses are **0.05 in the
contrast mean and 0.10 in every fold**. Near-zero target variance is unavailable.

Identity reproduces all moments/folds exactly. The library-PC-mean-only control
removes all within-library coordinate variation while retaining original RNA;
its within-library skill is zero. Baseline mean skill must exceed 0.05 to qualify
control sensitivity. Fourteen comparisons are sensitive and eighteen are not.
Even identity cannot qualify the complete gate when controls are insufficient.

Some baseline decoder R² values are negative. They mean this fitted decoder is
worse than each held-out library's own mean under this metric, not that the RNA
or PCs contain no signal. The fit mixes between- and within-library objectives;
decoder calibration requires a separate development experiment. These frozen
failures must remain intact if later models change that objective.

## Complete comparisons

H = 93 matched Hallmark members; F = the overlapping five-gene program.
Each method column is **mean loss / maximum fold loss**, baseline minus candidate.
Negative loss indicates improvement. A sensitive comparison fails when either
loss exceeds its margin. Insensitive rows remain measured but cannot qualify
biological preservation.

| Population | Treatment | Program | Baseline mean R² | Sensitive | Native 7 | Harmony 7 | Harmony 19 | Harmony 41 |
| --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: |
| Bcell | IFNa | H | -0.282533 | no | 0.008128 / 0.039082 | 0.006113 / 0.035715 | 0.006923 / 0.017548 | 0.006886 / 0.035578 |
| Bcell | IFNa | F | -0.609161 | no | 0.014352 / 0.048596 | 0.013209 / 0.047539 | 0.008424 / 0.035801 | 0.013834 / 0.045809 |
| Bcell | IFNb | H | -0.292730 | no | 0.014531 / 0.037481 | 0.012574 / 0.034601 | 0.013811 / 0.024245 | 0.013340 / 0.034189 |
| Bcell | IFNb | F | -0.613443 | no | 0.021392 / 0.064822 | 0.020112 / 0.064112 | 0.015611 / 0.053546 | 0.020956 / 0.062809 |
| Bcell | IFNg | H | 0.063957 | yes | 0.001141 / 0.015495 | -0.000762 / 0.013161 | 0.000245 / 0.010856 | -0.000431 / 0.013951 |
| Bcell | IFNg | F | 0.045481 | no | -0.000176 / 0.003724 | -0.000336 / 0.003766 | -0.001447 / 0.002652 | -0.000238 / 0.004040 |
| Bcell | IFN-L1 | H | 0.132824 | yes | -0.002616 / 0.016064 | -0.002975 / 0.015608 | -0.004028 / 0.011719 | -0.003009 / 0.015404 |
| Bcell | IFN-L1 | F | -0.005252 | no | -0.002985 / 0.005322 | -0.002671 / 0.005846 | -0.006631 / 0.002512 | -0.002477 / 0.005199 |
| Monocyte | IFNa | H | 0.461946 | yes | -0.005280 / 0.028542 | -0.007032 / 0.022505 | -0.005247 / 0.029638 | -0.005006 / 0.026884 |
| Monocyte | IFNa | F | 0.355779 | yes | 0.011814 / 0.046387 | 0.011068 / 0.039410 | 0.012040 / 0.048777 | 0.011523 / 0.044411 |
| Monocyte | IFNb | H | 0.463483 | yes | -0.005787 / 0.029962 | -0.007427 / 0.023190 | -0.005777 / 0.031376 | -0.005704 / 0.027870 |
| Monocyte | IFNb | F | 0.345129 | yes | 0.006927 / 0.049061 | 0.005948 / 0.040390 | 0.006759 / 0.051372 | 0.006206 / 0.046121 |
| Monocyte | IFNg | H | 0.415722 | yes | **0.053781 / 0.101709** | **0.055151 / 0.097165** | **0.053859 / 0.098655** | **0.054844 / 0.100033** |
| Monocyte | IFNg | F | 0.326662 | yes | **0.065771 / 0.116470** | **0.068350 / 0.113292** | **0.064036 / 0.116071** | **0.065390 / 0.116709** |
| Monocyte | IFN-L1 | H | 0.537780 | yes | -0.015149 / 0.046448 | -0.014726 / 0.050212 | -0.014504 / 0.040088 | -0.014891 / 0.045717 |
| Monocyte | IFN-L1 | F | 0.455269 | yes | **0.031582 / 0.252738** | **0.033485 / 0.263223** | **0.031081 / 0.236142** | **0.031794 / 0.252282** |
| NK | IFNa | H | -0.441691 | no | 0.033112 / 0.083439 | 0.032333 / 0.080794 | 0.030846 / 0.080300 | 0.032419 / 0.083035 |
| NK | IFNa | F | -0.959054 | no | **0.086238 / 0.148303** | **0.085848 / 0.147520** | **0.083741 / 0.145583** | **0.085452 / 0.154636** |
| NK | IFNb | H | -0.448468 | no | 0.032544 / 0.080058 | 0.031655 / 0.077333 | 0.030140 / 0.076444 | 0.031836 / 0.079377 |
| NK | IFNb | F | -0.983298 | no | **0.081948 / 0.150004** | **0.081387 / 0.148395** | **0.079173 / 0.147943** | **0.081103 / 0.157853** |
| NK | IFNg | H | 0.129715 | yes | 0.012328 / 0.027064 | 0.012873 / 0.029481 | 0.012438 / 0.027641 | 0.012464 / 0.027993 |
| NK | IFNg | F | 0.028259 | no | 0.000436 / 0.001955 | 0.000517 / 0.001951 | 0.000439 / 0.001926 | 0.000483 / 0.002044 |
| NK | IFN-L1 | H | 0.104062 | yes | 0.006270 / 0.028122 | 0.006651 / 0.030163 | 0.006448 / 0.029348 | 0.006500 / 0.029069 |
| NK | IFN-L1 | F | 0.028581 | no | 0.001227 / 0.006974 | 0.001250 / 0.006718 | 0.001130 / 0.006429 | 0.001223 / 0.006911 |
| Tcell | IFNa | H | -0.181607 | no | -0.010093 / 0.033627 | -0.010347 / 0.044465 | -0.009525 / 0.030134 | -0.010339 / 0.033703 |
| Tcell | IFNa | F | -0.560828 | no | -0.001648 / 0.079860 | -0.001690 / 0.097504 | -0.001736 / 0.075302 | -0.002040 / 0.082833 |
| Tcell | IFNb | H | -0.216580 | no | -0.006390 / 0.034761 | -0.006496 / 0.046938 | -0.005563 / 0.031428 | -0.006517 / 0.034724 |
| Tcell | IFNb | F | -0.638252 | no | -0.002803 / 0.071319 | -0.002860 / 0.088625 | -0.002499 / 0.067560 | -0.003148 / 0.073432 |
| Tcell | IFNg | H | 0.283946 | yes | 0.001921 / 0.004350 | 0.002171 / 0.004835 | 0.002045 / 0.003644 | 0.002140 / 0.004644 |
| Tcell | IFNg | F | 0.037452 | no | -0.000108 / 0.000578 | -0.000018 / 0.000726 | -0.000040 / 0.000851 | -0.000084 / 0.000725 |
| Tcell | IFN-L1 | H | 0.163441 | yes | 0.003554 / 0.006249 | 0.004079 / 0.006323 | 0.003409 / 0.005978 | 0.003845 / 0.006292 |
| Tcell | IFN-L1 | F | 0.042632 | no | 0.001907 / 0.008077 | 0.001924 / 0.008312 | 0.001760 / 0.007381 | 0.001904 / 0.008119 |

The three sensitive failures are Monocyte IFNg for both programs and Monocyte
IFN-L1 for the five-gene program. Each candidate passes both margins in 11/14
sensitive comparisons. Two additional raw margin failures occur in insensitive
NK IFNa/IFNb five-gene comparisons; they are retained without promoting them to
strong biological evidence. These observational proxy/engineering gates do not
by themselves establish causal biological damage.

## Numerical verification and runtime scope

Six focused tests pass, including exact mapping/coverage rejection, signed
sparse scoring, empty/constant targets, shuffled library moments, independent
weighted scikit-learn fitting, identity and erasure controls. Seven full data
evaluations complete in roughly ten seconds each with at most 132,300,800 bytes
peak RSS. Five independent checks then fit every actual fold/program directly
from weighted per-cell rows using augmented SVD least squares, rather than
library sufficient statistics. They compare coefficients and direct held-out
residuals for baseline, native and all three Harmony matrices. All pass; maximum
R² discrepancy is 7.15e-13 and maximum coefficient discrepancy is 3.10e-13.
Each check takes about 12.3 seconds, with at most 683,720,704 bytes peak RSS.
These are concurrent-host observations, not a controlled speed comparison.

The native matrix is the already published [seed-7 integration](FULL_INTEGRATION_RESULTS.md),
executable SHA256
`e360645362ac0707cb497eaa64277e5bf792b28d79bddea707421fdcd043495d`.
It predates sequential-ridge optimization. All three [Harmony matrices](HARMONY_REFERENCE.md)
are the original complete seeds 7/19/41. No integration refit is claimed here.

This is an independent measured-program reference, not a newly executed native
million-cell program publication. The existing native pseudobulk JSON report
is 526,518,683 bytes against its 536,870,912-byte limit, leaving insufficient
space for the added per-cell program arrays. Native full-scale scoring needs a
compact separate program artifact; increasing the global report limit does not
resolve that ownership/storage gap. Existing smaller native program evidence
remains separate. Rare-cell, neighborhood, annotation and unseen-context gates
also remain open.

## Archive and reproduce

The [archive manifest](evidence/2026-09-10-integration-programs/manifest.json)
retains all per-cell program scores, detections and totals; exact definitions,
protocol, design, execution sources/freezes; seven results with every fold and
library moment; five dense-oracle results; logs, test outcomes and original
candidate receipts. Original H5AD, metadata and five latent matrices remain
external with exact byte counts/hashes and the existing study locations in
recorded commands. Historical QC receipt status is preserved unchanged.

```sh
python Tools/Omics/Benchmarks/HIRISA/verify_archive.py \
  Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-10-integration-programs
python -m unittest discover -s Tools/Omics/Benchmarks/HIRISA \
  -p test_integration_programs.py
```

`prepare_program_reference.py`, `check_integration_programs.py` and
`check_program_dense_oracle.py` expose their required original inputs through
`--help`. Frozen coordinator scripts in the archive retain every executed
command. `archive_integration_programs.py` preserves completed evidence without
claiming external matrix re-execution.
