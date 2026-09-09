# Untreated-cell correlation descriptors: prospective second experiment

Specified 2026-09-09 after the co-response experiment, before scoring this model.
The first experiment informed this design. This is a reused research benchmark,
not a new untouched external validation set; do not tune this protocol after scores.

Keep the same 105 held-target single-condition folds and full feature universe.
All conditions involving a held target remain excluded from fitting. The only
change to target information is a descriptor derived exclusively from the 11,855
untreated cells. No perturbed-cell count enters descriptor selection or estimation.
Use the byte-pinned H5AD and verified reference cell library totals.

Stream every source gene, selecting control rows before normalization. Transform
nonzero counts to log1p(CPM) using each control cell's full library. Rank eligible
landmarks by population log-expression variance (Ensembl ID ties), requiring
expression in at least 10 control cells. Exclude all exactly matched intervention
target genes from landmarks. Keep the first 2,000. This is a declared variance
selector, not Scanpy's HVG flavor or a fitted regulatory network.

A target descriptor is its Pearson correlation with each landmark across control
cells. Targets require a unique exact source-symbol match, positive control
variance and expression in at least 10 controls. Missing/weak targets remain
explicitly unsupported, without an expression-response fallback labeled as a
specific prediction. Every target is still scored for no-change and mean-single.
Compute moments and cross-products with sparse control matrices; never construct
a dense cells-by-all-genes matrix. Independently check selected correlations.

For each held target, use descriptors of other supported targets to predict
whole-transcriptome single-target log1p-CPM responses by ridge, intercept enabled,
alpha=1. Center and population-standardize descriptor coordinates on training
targets only; constants use 1; divide by sqrt(2,000). No held-response tuning.
Compare dual solve with independent SVD Ridge. A fixed one-step rotation of
sorted training responses supplies the matched shuffled control, not a p-value.

Predictors: noChange, meanSingleResponse (all other 104 singles),
controlCorrelationRidge, shuffledControlCorrelationRidge. Clip negative predicted
log-expression, retain implied CPM sums, no reclosure. Freeze predictions before
separate scoring. Use all genes and each fold's training top-1,000 expression
panel; retain all failures and compare methods only on matched supported targets.
Report correlation/sign/MAE as well as RMSE. Repeat the complete descriptor and
prediction preparations, and verify source-control identities, normalization
sums, fold exclusions and replays. Hash protocol, source, helper and outputs.

Control coexpression is observational and can reflect cell state or technical
confounding. It does not establish causal influence or a mechanistic connection.
No donor, unseen-pair, uncertainty, native-runtime or tissue qualification follows.
