# Frozen mouse sample annotation comparison, 2026-09-21

Use only observational cell-type annotation. Training is the complete original
GSE84133 GSM2230761 mouse1 count table; query is the complete GSM2230762 mouse2
table. Preserve all source cells and genes. Source author annotations are the
scoring target, not independently re-established biological identities. The
source contains two mouse strains; this is one sample/strain holdout and is not
independent human-donor, laboratory, disease, intervention or phenotype validation.
The four human Baron donors informed earlier software development. No prior
mouse classifier result was found in the repository, but public labels existed
before this run. Blinding is limited to this execution, not prospective biology.

Freeze source/member SHA-256, exact runtime, plans and this protocol before fit.
Write training labels only into training H5AD. Query H5AD contains no labels;
its preparation discards the label field without retaining or inspecting it.
The reveal reads the original mouse2 labels only after native predictions and
all baseline predictions are written and hashed. No post-reveal model changes.

Native NumiVivo balanced multinomial logistic regression uses the existing
settings: library target 10,000, 2,000 training-selected HVGs, 20 mean bins,
20 centered PCs, maximum Krylov basis 128, PCA residual tolerance 1e-6, seed 7,
penalty 1, gradient tolerance 1e-8, at most 2,000 iterations and 20 billion
classifier work units. Scaling, feature selection and PCA use training only.
Each query normalizes by its own complete RNA library; no query cohort fit,
batch correction, integration or label-derived feature selection is allowed.

Baselines are training-majority (sorted-label tie break), nearest class centroid
in the same frozen unstandardized PCA coordinates (squared Euclidean distance,
sorted-label tie break), and native uniform kNN with k=15. No baseline sees test
labels. All methods predict every query cell. Empty query libraries receive an
explicit unmapped label and count as errors; missing genes or invalid counts
abort rather than being rounded, zero-filled or omitted.

Primary metric is macro-F1 over the union of training classes, revealed query
classes and predicted classes, with zero division equal to zero. Unseen query
classes remain in the metric and confusion matrix with zero recall if missed;
no class merging, rare-class removal or favorable subset is permitted. Success
requires native logistic macro-F1 at least 0.02 higher than BOTH centroid and
native kNN, and strictly higher than majority. Each cell contributes once;
cell counts are not biological replicate counts. Accuracy, balanced accuracy,
per-class support/precision/recall/F1 and confusion matrices are secondary.
No multiple-endpoint substitution or confidence claim from one held-out sample.

Native reconstruction and a repeated mapping must reproduce results. Independently
recompute softmax and training centroids from retained native parameters/scores.
Retain failures, exact wall time, runtime hash and actual CPU backend. This
experiment does not qualify Metal acceleration or the current complete package.

Primary source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133
and the linked GSM2230761/GSM2230762 source tables. Original data: Baron et al.,
Cell Systems (2016), DOI 10.1016/j.cels.2016.08.011.
