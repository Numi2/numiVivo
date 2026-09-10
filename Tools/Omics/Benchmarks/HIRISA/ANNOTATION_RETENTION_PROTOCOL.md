# HIRISA donor-held-out annotation retention protocol

Declared 2026-09-11 after the existing HIRISA response, program, clustering and
preparation-transfer outcomes were inspected, before this annotation-retention
experiment. This is a diagnostic on an inspected study, not independent biological
validation or a prospective cell-identity predictor. Original author `celltype.l2`
predictions are comparison labels, not ground truth. No label-confidence filter,
expression-based inclusion, relabelling or integration refit is permitted.

Use all 1,612,594 original cells, the exact 131-library metadata, and the existing
20-dimensional baseline, native seed-7 and Harmony seed-7/19/41 matrices. Validate
original row/barcode/library/donor order. Retain every original preparation ×
treatment stratum. Hold out each donor within each stratum, pooling its technical
libraries without treating them as independent donors. Every cell is a query
exactly once. Record unavailable training, unseen classes, all confusion entries
and every class with positive query support; do not omit unsuccessful strata.

Fit fixed class-balanced multi-output ridge to one-hot training labels. Each
training class has equal total weight; cells within that class have equal weight.
Use training-only weighted centers and standard deviations, all twenty PCs,
fixed ridge 1 and an unpenalized intercept. No query labels enter fitting.
Prediction selects maximal fitted score; exact ties split confusion weight
uniformly among tied classes. No class probabilities or calibrated uncertainty
are claimed. Classes absent from training receive no predictions and zero recall
when present in the query. An all-zero-PC erasure control must produce uniform
votes over training-supported classes. Exact identity must reproduce results.

The primary retention endpoint is each stratum/class's equally weighted mean
held-out-donor recall, with maximum individual-donor recall loss also reported.
A class is operationally rare when its complete metadata frequency is below 1%
within its preparation/treatment stratum; this designation uses no RNA or score.
All classes are evaluated, with rare classes reported separately. Query/train
support below twenty cells is measured but insufficient for qualification. A
comparison requires every positive-query donor fold to have sufficient support
and baseline recall at least 0.05 above the erasure control in the donor mean.
Allowed losses are at most 0.05 in mean recall and 0.10 in every fold. Missing or
insensitive comparisons cannot qualify complete preservation. These declared
engineering margins do not by themselves prove biological damage or preservation.
Report precision/F1, full confusion, overall accuracy and balanced recall as
secondary descriptions; no averaging over only successful comparisons.

Compute statistics in bounded row blocks with no dense cells-by-genes array.
Reconstruct every actual training fit independently from class-weighted per-cell
rows with augmented SVD least squares, then compare query confusion. Verify each
matrix's original payload hash and explicit coordinates. Bind source, metadata,
protocol, scripts, label ledger and all results with hashes. Preserve prior
failures and the unresolved independent-study/rare-cell/annotation qualifications.
