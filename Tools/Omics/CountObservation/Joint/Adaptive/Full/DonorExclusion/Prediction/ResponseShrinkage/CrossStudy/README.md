# Native response shrinkage: cross-study transfer remains unqualified

The training-selected response-shrinkage rule from the 13-fold development
experiment now has an executed Swift implementation. All 26 matched-panel
within-study and cross-study fits/predictions complete. Independent NumPy
calculations reproduce every penalty's inner-donor validation loss, every final
coefficient and all 308,984 predicted values (maximum prediction discrepancy
3.67e-15). Python prepares normalized donor pseudobulk inputs; this qualification
covers native fitting and prediction, not an end-to-end native H5AD command or a
full product build.

## Complete result

The panel is the historical exact, unique symbol correspondence: 11,884 genes.
Normalization retains every source feature in each library denominator. Scores
are equal-donor mean RMSE on log1p treated CPM, equivalently response RMSE with
observed query control. These values differ from the prior restricted NB-eligible
within-study panel and its pooled gene RMSE; do not merge their denominators.

| Training → query | Mode | Candidate RMSE | Training mean | No change | Gain over mean |
| --- | --- | ---: | ---: | ---: | ---: |
| HIRISA → Kang, 8 donors | Transfer | 1.130197 | 1.135572 | 1.220670 | 0.47% |
| Kang → HIRISA, 5 donors | Transfer | 1.651690 | 0.521068 | 0.341591 | -216.98% |
| Kang → Kang, 8 omissions | Within study | 0.941749 | 1.156071 | 1.220670 | 18.54% |
| HIRISA → HIRISA, 5 omissions | Within study | 0.104956 | 0.106039 | 0.341591 | 1.02% |

HIRISA → Kang passes the historical strict-better-than-both aggregate rule but
misses the separately reported 5% practical gain threshold; one donor loses to
training mean. Kang → HIRISA fails both baselines and all five donors lose to
training mean. Within-study HIRISA loses one donor on this complete panel. This
is not a qualified transferable replacement. No dataset-specific winner mixture
is promoted. Original joint, transport and cross-study ridge results remain
unchanged. Further outcome-driven tuning on these cohorts cannot provide fresh
validation.

The studies differ jointly in donor health, treatment duration, assay and cell
preparation. Kang's annotated B-cell and HIRISA's enriched-B-cell preparations
are not equivalent populations. These comparisons do not identify a causal
reason for the transfer failure or qualify phenotype/clinical prediction.

## Method and provenance

For each gene, regress training log1p treated-minus-control CPM response on
training log1p control CPM. Shrink the centered slope with penalties
`0, .01, .1, 1, 10, 100, infinity`, where infinity is training mean response.
Choose one penalty per fold by inner donor mean squared error on all panel
genes; exact ties prefer stronger shrinkage. Refit on all admitted training
donors and clip the predicted treated log1p CPM at zero. Both full-study transfer
models select 0.1. No query-treated outcome enters fitting or prediction.
Unlike the previous count-kernel experiment, this full-panel transfer test has
no NB availability filter. It is an observed-pseudobulk regression with no new
NB likelihood, parameter uncertainty integration or calibrated interval.

`protocol.json` freezes the equation, selection, endpoint and native/source
hashes before execution. All training/query H5AD and plan files must match the
historical input freeze. `prediction-freeze.json` binds all 26 input/output pairs
before scoring. The scorer independently reconstructs normalization from the
original source count arrays and checks every training and query input vector.
Repeated scoring is byte-identical. This uses already inspected studies and
transfer outcomes; it is development evidence, not an untouched external test.

`Tests.swift` passes four analytical cases (mean limit, known slope, clipping,
constant control) and seven identity/dimension/penalty/numeric rejections.
`build-attempt1.log` retains an initial Swift expression type-check failure;
splitting that expression resolved compilation before execution.

## Reproduction

`ResponseShrinkage.swift` is a scoped experimental implementation, retained here
with its driver rather than promoted to a production model interface. Compile:

    swiftc -O ResponseShrinkage.swift Driver.swift -o response-shrinkage
    swiftc -O ResponseShrinkage.swift Tests.swift -o tests
    ./tests

Use the scientific Python environment with AnnData and NumPy. Restore the
historical `numivivo-cross-study-ifnb-20260911/inputs` directory and adapt the
explicit paths if necessary. Run `run.py`, `verify.py`, then `score.py`. Native
output paths must not already exist. The archive retains exact native sources,
binaries, all normalized model inputs and outputs, protocol, hash bindings,
verification, tests, scores and logs; historical H5AD/count arrays remain in their
original study. It does not duplicate the atlas or claim a Metal speedup.

The [observed-range and error decomposition](Support/README.md) shows that Kang → HIRISA fails even inside the training range in every donor (1.82–1.88× mean-baseline RMSE). Greater range coverage does not identify the more accurate direction. A range filter is not a validated fix; the next selection design needs held-out study/context structure, not further tuning on these query outcomes.
