# Three-study held-out response model: transfer still fails

A native Swift model now selects shrinkage using held-out training studies,
then predicts a third study from query controls only. All 75 donor predictions
complete on the existing 11,800-gene shared panel. An independent NumPy
implementation reproduces every study weight, every penalty's validation loss,
all model coefficients and all 885,000 predicted values. Separate scalar scoring
checks all 225 donor/baseline RMSE values.

## Results

Each row reports equal-donor mean RMSE on log1p treated CPM. The mean baseline
uses equal weight per training study, then equal weight per donor within study.

| Entire study held out | Donors | Selected penalty | Candidate | Balanced training mean | No change | Gain over mean | Gain over no change |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Kang | 8 | 10 | 1.120253 | 1.135852 | 1.224197 | 1.37% | 8.49% |
| HIRISA | 5 | 1 | 0.684251 | 0.336026 | 0.342458 | -103.63% | -99.81% |
| GSE181897 | 62 | 100 | 1.067425 | 1.069742 | 1.121616 | 0.22% | 4.83% |

All Kang and GSE181897 donors improve against both baselines. All five HIRISA
donors lose against both. No study reaches the declared 5% gain over both
baselines; the model is not promoted. The HIRISA failure persists despite
stronger shrinkage chosen without HIRISA training outcomes. Do not pool the
three rows to hide that failure or select a different model by query-study
outcome. Prior two-study values use a different panel/training contract; do not
attribute their numerical differences solely to study-level selection.

This is development reuse of three previously inspected cohorts. In particular,
GSE181897 now enters training in two outer folds, explicitly superseding its
query-only role for this new experiment. Its original external test protocol,
results and failures remain unchanged. No untouched external, prospective,
clinical or calibrated-uncertainty claim is made. Source study exclusion is
verified; absence of participant overlap across source studies is not established.

## Training and selection

For each gene use x = log1p(control CPM), d = log1p(treated CPM) - x.
Within any training set, each study has equal total weight and each of its
donors has equal weight. With weights summing to one:

    slope = weightedCovariance(x, d) /
            (weightedVariance(x) + lambda * (1 - sum(weights²)))
    prediction = max(0, queryControl + weightedMean(d)
                        + slope * (queryControl - weightedMean(x)))

A zero denominator gives zero slope; infinite penalty is the mean-response
limit. The adjustment recovers the previous `(n-1)` ridge denominator when all
n donors have equal weight. Penalties are fixed at 0, .01, .1, 1, 10, 100 and
infinity. For each outer held-out study, leave each of the two training studies
out in turn, compute all-gene/equal-donor MSE within the validation study, then
average the two study losses equally. Choose one penalty for the complete
outer fold; exact ties prefer stronger shrinkage. Refit on both training studies.

All original source RNA features remain in normalization denominators. The
panel uses the already qualified exact symbol correspondence. No alias guessing,
zero padding, gene filtering by outcome or new cell-label fitting occurs.
GSE181897's primary author code establishes B=IFNB and C=control. All 62 paired
donors are retained; the original two donors missing IFNB remain excluded.
Source count aggregation is reconstructed again against every native GSE181897
pseudobulk row before preparing this experiment. HIRISA is enriched B cells,
whereas the other endpoints use author B-cell annotations; assay, treatment
time/dose, population and study remain confounded.

## Verification and reproduction

`StudyShrinkage.swift` is a scoped experimental fitter/predictor. Python prepares
normalized donor pseudobulk inputs; this is not a new NB observation likelihood,
full-product CLI integration, end-to-end native H5AD qualification or Metal
speedup. Seven software tests cover unequal donor counts with equal study
weight and invalid study/donor/shape/numeric inputs. The real study is the
biological evidence; those edge cases do not qualify biology.

`protocol.json` binds method and source before execution. `input-freeze.json`
binds original source artifacts, normalized cohorts and all three native input
files. Query-treated values are absent from each native input; only the scorer
uses the corresponding frozen query outcomes. `prediction-freeze.json` precedes
scoring, and `pre-score-bindings.json` retains the scripts and verification used.
Independent maximum prediction discrepancy is below 3.56e-15.

Compile with `swiftc -O -parse-as-library StudyShrinkage.swift -o study-shrinkage`.
Use the scientific Python environment with AnnData, SciPy and NumPy. Run
`prepare.py`, `run.py`, `verify.py`, `test.py`, `score.py`, then `verify_scores.py`.
Adapt explicit local paths when restoring; new output directories must not
already exist. The evidence archive retains sources, binary, normalized
cohorts, native inputs/outputs, protocols, checks and scores. Original source
H5AD/count files remain in their retained study directories and are not duplicated.

The next scientific gap is how to account for differences between study contexts
in the response model, rather than further donor-only penalty tuning. This test
establishes that the current study-balanced linear response model is insufficient;
it does not identify which assay, exposure or population difference causes failure.
