# Frozen held-out-animal count prediction

Use the complete author-deposited Crowell post-QC count scope already verified
against Scanpy: eight independent mice, four Vehicle and four LPS. Evaluate all
16 choices of one held-out animal per arm in all seven replication-eligible
populations. CPE remains unavailable. Training has three animals per arm; no
cells or genes are selected using held-out effects or significance.

Fit the existing native NB cohort owner on training pseudobulks only, using the
existing count/presence/support gates, Gamma parametric dispersion trend and
library-size normalization. The normalization change from the earlier
median-ratio comparison is fixed before any held-out score. All normalization,
dispersion, gene availability and empirical-prior estimation use training
counts. Compare unshrunk MLE, fixed normal contrast prior SD=1 log2 unit and
empirical weighted-quantile normal contrast prior. All three share the exact
training final dispersions; nuisance coefficients are jointly refitted.

Fit input files contain only the six training pseudobulk rows. The independent
scoring command receives the frozen fitted model and two held-out rows. Test
offset is log(observed test library total) minus the mean log training library
total. Thus the score conditions on known library depth; it is not prediction
of absolute sequencing yield. This is plug-in NB predictive scoring, not an
integrated posterior predictive distribution or calibrated interval assessment.

Primary comparison: empirical versus unshrunk mean held-out NB log score,
averaged equally over genes admitted by training, then over the four predictions
for each actual animal, then equally over eight animals and seven populations.
Positive delta is descriptive improvement. Fixed-versus-unshrunk and
empirical-versus-fixed deltas and squared log1p normalized-count error are
secondary. Retain every animal/population result, failed fit and unavailable
gene. If either MAP is unavailable for a training-tested gene, that fold's
primary comparison is unavailable; do not silently narrow its gene family.

Report gene coverage and training-only filter/status counts in each fold.
Folds overlap, as do animal measurements across populations; do not treat 112
folds, individual genes or populations as independent animals, and do not attach
an independence-based p-value or confidence interval. This previously inspected
study cannot establish prospective cross-study benefit, FDR or posterior
coverage. Do not tune any method or change this protocol after viewing scores.

Execution uses a bounded native measurement harness compiled from the unchanged
production Omics source files. It calls the existing cohort and MAP owners; it
does not implement a second estimator. First compare its full-cohort output
with the previously qualified product report, then run training-only folds.
Independently verify all conditional scores, information and held-out log masses
with NumPy/SciPy. Keep the existing product CLI, source aggregate and executable
hashes distinct from this measurement harness and its compact output format.
