# Native balanced multinomial reference mapping protocol

Declared 2026-09-11 before native logistic fits or new outcome scoring. The
previous complete Baron results have been inspected: balanced logistic improves
macro-F1 over uniform kNN in all four folds, with retained rare-class failures.
This is implementation and full-cohort requalification on inspected development
data, not an untouched biological validation or model-selection experiment.

Use all four original donor-held-out Baron splits, every source gene and every
query cell. Retain exact frozen source/label identities and the original sparse
HVG/PCA settings. Query labels are unavailable to inference. Fit training-only
score standardization and class-balanced multinomial logistic regression with
L2 coefficient penalty 1 (equivalent to C=1), unpenalized intercepts, and weights
n/(number_of_classes * class_count). Do not tune the penalty, classes, donors,
features or stopping criterion using held-out labels. Preserve uniform kNN as
the default and as the existing comparison.

Native optimization must report its objective, gradient convergence and charged
work. Nonconvergence is an error, not a qualified model. Independently check
standardization, class weights, loss/gradient, all query probabilities and labels
against scikit-learn/SciPy on identical native training coordinates. Compare the
original complete benchmark predictions and per-label results, retaining any
solver differences. The numerical target is maximum probability difference
1e-4 against a tighter independently fitted optimum, with label disagreements
and margins explicitly reported; numerical agreement is separate from accuracy.

Freeze native predictions before scoring author labels. Report every donor's
accuracy, macro-F1, per-class precision/recall and confusion matrix, including
zero-recall classes and worse-than-kNN results. Do not claim calibrated class
probabilities, novel-class rejection, authoritative cell identities or independent
study/tissue transfer. Single-cell counts remain sparse; dense state is limited
to PCA scores, model parameters and class outputs. Verify source reconstruction,
exact replay, default kNN compatibility, empty queries, label-leak/overlap and
work-budget rejection, and tampered models.

Execute complete native fits/maps on the physical Mac mini after checking active
workloads. Retain all source/plan/executable hashes, failures and complete
restorable bundles. Reuse existing frozen donor inputs and verified completed
work; bound storage by compressing finished run-owned scratch only after archive
byte and open-handle checks. No original data or previous evidence is disposable.
