# Baron donor-held-out reference mapping protocol

Predeclared before running the benchmark, 2026-09-09.

Use all original Baron human UMI cells and features. Evaluate four leave-one-donor-out folds. Author cell types are reference labels; held-out labels are evaluation-only. Keep all rare labels and report absent training classes. Human4 disease is confounded with donor: this does not evaluate disease effects.

Native H5AD projection creates training and query bundles retaining the full gene axis. Native streamed reduction fits only training cells: library normalization target 10,000, Seurat dispersion HVG 2,000, 20 mean bins, centered PCA 20 components, Krylov basis 128, residual tolerance 1e-6, seed 7. Scanpy independently verifies training HVGs and PCA eigenvalues. Query normalization uses the complete original feature universe, then projects selected log counts using training means and loadings. No query fitting or integration.

Two fixed external reference classifiers: uniform 15-neighbor Euclidean voting (brute force), and balanced multinomial logistic regression (C=1, lbfgs, max_iter=2000, tol=1e-6) on training-standardized PC scores (population standard deviation). No hyperparameter selection on query labels. Predictions are candidates, not authoritative labels; probabilities are uncalibrated.

Retain source/binary hashes, source axes, fold membership, training centers/loadings, classifier parameters, predictions, per-label precision/recall/F1/support, confusion matrices, accuracy, balanced accuracy and macro F1. Preserve warnings and failures. These results establish an external target for a subsequent native frozen-reference implementation; they do not qualify native learned annotation or general cross-study mapping.
