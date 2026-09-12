# Frozen response-weight audit

The failed Parse context model already balances its three training studies:
GSE181897 has 62 donors, HIRISA 5 and Kang 8, but each receives one third of base
weight. Donor count imbalance is therefore not an explanation for base weighting.

All twelve fitted query response-weight vectors are nonnegative and sum to one
within 1e-12. Negative weight mass is zero. Constraining negative weights to zero
and renormalizing would not alter these predictions. This rules out that particular
correction for this frozen model; it does not establish accurate response transfer.

| Study | Base weight | Fitted total weight across the twelve queries |
| --- | ---: | ---: |
| GSE181897 | 1/3 | 0.4021–0.4080 |
| HIRISA | 1/3 | 0.0495–0.1087 |
| Kang | 1/3 | 0.4846–0.5442 |

Nearest training-control donors by full-panel RMS distance are Kang patient_1015
or patient_1039. Control similarity does not imply that treated responses transfer.
The model's inspected error remains excessive response magnitude plus imperfect
direction. These observations support investigating exposure/population-dependent
response mapping; they do not identify a proven correction or validate a new model.

`run.py` reopens the hash-verified frozen input and native output, checks equal study
base weight and weight sums, and records every query in `audit.json`. It reads no
Parse target file, fits nothing and changes no predictions. The overall investigation
is retrospective, motivated by inspected outcomes; it is not independent validation.
The source hashes match the prior failure diagnosis. Original files remain at
`/Users/n/numivivo-parse-context-evaluation-20260912`, with this audit at
`/Users/n/numivivo-context-weight-audit-20260912`. The driver requires NumPy and a
fresh output directory for reproduction.
