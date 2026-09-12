# Complete feature-level temporal error audit

The earlier-duration-only model passes its predefined donor-average gate, but
per-feature error is not uniformly improved. This retrospective diagnostic weights
donors equally and available times equally within donor before comparing each
feature's mean squared error across the original 18 observations.

| Baseline | Genes with lower candidate MSE | Higher MSE | Tied within absolute 1e-12 |
| --- | ---: | ---: | ---: |
| Training mean | 9,202 | 3,769 | 22 |
| No change | 9,662 | 3,308 | 23 |

All 12,993 features remain, including constant and weakly expressed features. The
absolute tie tolerance is numerical bookkeeping, not statistical significance or a
biological effect threshold. Gene outcomes are correlated; these counts are not
independent replications or an estimate of general gene-level accuracy.

The complete compressed CSV retains all three errors for every feature. `audit.json`
records the largest regressions, source/frozen-prediction hashes and table hash;
`audit.py` verifies every native output hash before calculation. A separate CSV
check reproduced the counts and verified all unique feature identities. No features
were selected for a new model, no output changed, and no acceptance rule was revised.
This does not supersede the donor-level result; it bounds its interpretation.
Additional donors, independent studies and calibrated gene-level uncertainty remain
necessary. Runtime files are at `/Users/n/numivivo-duration-feature-audit-20260912`.
