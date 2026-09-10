# Independent profile and reference-boundary audit

Declared after the dispersion-stage diagnosis at 0c14113 and an exploratory
Hagai seed-1 gene check found a higher interior objective than the reported
DESeq2 boundary. This is post-score diagnosis, not independent calibration.

For every one of the twenty original null splits, reproduce DESeq2's original
gene-wise estimates and retain its actual pre-dispersion fitted means, iteration
counts and bounded initial estimates. Evaluate its own grid dispersion routine
on all eligible genes with exactly those means, counts, design and Cox–Reid
correction. Preserve original estimates and grid alternatives separately.
Do not replace the frozen comparison, refit its downstream trend, select a new
BH family or claim that a higher likelihood establishes better calibration.

Independently check the likelihood difference between original and grid
estimates for every eligible gene. Use the exact integer-count rising-factorial
identity for NB log mass; it avoids subtracting huge log-gamma values near the
Poisson boundary. Use QR for the information determinant. Count improvements
above 1e-4 log-likelihood units, while retaining every signed difference.

For native numerical checks, take the sixteen genes with largest native
gene-wise dispersion among native-interior/DESeq2-below-1e-6 disagreements in
seed 1 of each study, plus all twenty original Hagai sham calls. Deduplicate
study/seed/gene identities. Record selection before profile evaluation. Compile
the exact native QR/NB owner sources and inspect native-selected, reference-
selected and grid-selected points, plus 41 logarithmically spaced dispersions
over the native bounds [1e-8,100]. Retain per-point failures. Independently refit
coefficients with statsmodels and compare stable likelihoods, means, effects,
standard errors and adjusted profiles. Also compare the native selected profile
against its surrounding independent grid, without claiming global optimality
for arbitrary inputs.

No native inference option or production default changes based solely on this
audit. Any actual numerical defect requires a source repair and new qualification.
