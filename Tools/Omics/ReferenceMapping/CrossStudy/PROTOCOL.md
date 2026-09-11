# Complete Kang–Ding reference transfer

Frozen before fitting this reference-transfer experiment. Both cohorts have
already been inspected for integration; Baron was used for classifier development.
This is a development transfer test between separately collected studies, not
an untouched prospective study or independent biological identity experiment.

Run both directions with the existing prepared original-count sources. Kang
training retains all 24,673 labelled cells, both conditions and eight donors.
Ding training retains every assigned source label (29,411 cells); 46 Unassigned
and 14,574 unavailable labels cannot supervise a classifier. Both queries retain
every cell: Ding 44,031 and Kang 24,673. No outcome-driven cell/gene selection.
Keep the original full gene axes and all RNA counts for library normalization.

The feature panel is the 14,976 exact Kang symbols appearing once in the original
Ding symbol column. Nineteen ambiguous Kang symbols and 711 absent symbols cannot
enter the panel. For Ding, only those unambiguous features receive their literal
source symbol as an explicit native feature ID; other IDs stay original. Preserve
the complete reversible ID table and original source. Do not resolve synonyms,
merge duplicates or interpret unmeasured genes as zeros. Every panel gene must
be measured by training and query. No cell labels or counts select the panel.

Train-only HVG binning over the panel requests 2,000 genes, 20 mean bins, 20 PCs,
128 basis vectors, residual tolerance 1e-6 and seed 7. Normalize each complete
source library to 10,000, then log1p. Fit native training-only standardization and
class-balanced multinomial logistic regression with penalty 1, gradient tolerance
1e-8, maximum 2,000 iterations and 100 billion charged objective/gradient terms.
No integration, hyperparameter search, probability calibration or query fitting.
Preserve every failed attempt and do not weaken numerical gates after inspection.

Freeze native models/predictions and reconstruct them before scoring. Independently
check full source axes/count totals, panel-only Scanpy HVG selection, PCA covariance
residuals, original-library query projection, standardization, objective, gradient
and every probability. Require probability agreement with independent tighter
scikit-learn fitting within 1e-4, explicit softmax within 1e-12, and normalized PCA
residuals within 1e-6 plus 1e-8 numerical comparison allowance. PCA signs are not
biological discrepancies. Compare a training-majority baseline on identical queries.

Keep the original eight/nine label vocabularies in native models and all rectangular
source-label versus prediction tables. They are not identical taxonomies: notably
Cytotoxic T cell is not synonymous with CD8 T cells; pDC is absent from Kang's
separate labels. Never invent a fine-label accuracy from these correspondences.
Report explicitly coarsened six-family diagnostics (T, B, NK, monocyte, dendritic,
megakaryocyte), grouping each study's stated subclasses only for evaluation.
Report family accuracy, balanced accuracy, macro-F1, per-family precision/recall,
confusion, label coverage and each original library separately. The descriptive
family target is balanced accuracy and macro-F1 at least 0.8 in both directions,
with recall at least 0.5 for each supported family. This analyst-declared coarse
target cannot qualify fine labels, rare classes, novel-class rejection or tissue
biology. Unavailable labels remain unavailable, with their full predictions retained.

Ding's pbmc1/pbmc2 are biological samples with unreported donor IDs; methods are
technical libraries. Kang donors/conditions remain explicit. No cell/technical
library count is treated as independent donor replication or clinical evidence.
Results include both Kang conditions and all Ding libraries without selecting the
best subset. Report failures even if numerical reconstruction passes.
