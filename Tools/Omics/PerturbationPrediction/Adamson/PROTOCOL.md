# Adamson UPR screen: fixed target-kernel validation protocol

Declared 2026-09-10 before downloading, preparing or scoring Adamson responses.
The preceding Norman implementation and development results are frozen at
`6a83c08d0708d3236db7f3136ddd4ab92119b18b`. This validation must not change lambda,
annotation inclusion, feature panels or primary comparisons after outcomes are
scored. A failure remains a failure of this fixed protocol.

## Study and inclusion

Use the complete published scPerturb H5AD for
`AdamsonWeissman2016_GSM2406681_10X010`: the UPR screen, GEO GSM2406681/GSE90546,
[Adamson et al., Cell 2016](https://doi.org/10.1016/j.cell.2016.11.048),
PMID [27984733](https://pubmed.ncbi.nlm.nih.gov/27984733/).
The [primary GEO sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2406681)
identifies human K562 cells, Cell Ranger UMI-collapsed counts and hg19 alignment.
The source is a different CRISPRi experiment from Norman's CRISPRa screen.
Pin the release version, published checksum, size, independently computed SHA-256,
processing code and exact condition/guide identities before fitting.

Keep all source cells and features in the audited aggregate reference. Record
unassigned/ambiguous guides and combinatorial conditions explicitly. Train and
score every unambiguous single-gene condition with a positive count library;
controls are the source-declared negative controls. Multiple guides for the same
gene belong to one held target. Never split guides for that target across train
and test. Exclude every condition containing the held target and all multi-target
conditions from fitting. Exclude unresolved perturbation assignments with counts
and reasons, without selecting targets based on expression effect or score.
Pooled cell/technical batch data do not establish independent biological replication.
Do not silently combine distinct biological cell contexts.

## Frozen method and annotations

Fit within this study from control and the other single-target aggregates. This
is independent-study validation of the algorithm, not transfer of Norman-trained
response coefficients or evidence for unseen cell types/donors.
Use the existing native `VivoTargetKernel` method unchanged: log1p-CPM with full
source gene denominators, direct-term Jaccard kernel, ridge lambda=1 with an
unpenalized intercept, clipping at zero and no CPM reclosure. Retain no-change,
all-training-single mean, supported-training-single mean and the same one-position
supported-response shuffle. No inner or outer hyperparameter selection.

Use the same GO rule as the frozen Norman protocol: exact original Ensembl gene
identity, MyGene human taxon 9606 verification, distinct direct BP/MF/CC terms,
all evidence classes, exclude NOT qualifiers and GO roots, no ancestry expansion
or alias guessing. Reuse exact matching captured identities when provenance and
service build agree; capture any new responses and service metadata. Missing,
ambiguous or retired identities and no-data annotations remain unsupported.
Do not use measured post-perturbation RNA as target descriptors. Retain annotation
citation overlap with PMID 27984733 and describe that current knowledge need not
be temporally independent of the experiment.

## Separation, verification and acceptance

Freeze annotations, fold selections and predictions before the scoring step.
The native fit receives only control plus training-single counts and descriptors.
The query plan receives only context and independent target descriptors. Outcomes
are read by the separate scorer. Mutating excluded counts must leave selected
training counts unchanged. Verify native model/prediction reconstruction, original
feature identities/counts against independent sparse HDF5/SciPy aggregation, all
native vectors against an independent centered-kernel scikit-learn fit, and exact
repeated scientific output. No dense cell-by-gene array.

Primary endpoint: mean per-target all-gene log-response RMSE on the same supported
held targets, comparing fixed GO with the all-single mean, supported-single mean,
and matched fixed shuffle. A useful GO-specific gain requires lower primary RMSE
than all three comparators. Report every target and worse-than-no-change failures.
Keep all available generic-baseline targets as a separate coverage cohort. Also
report the same training-selected top-1,000-mean-CPM panel, MAE, correlation and
sign metrics; this secondary panel cannot rescue a failed primary endpoint.

This protocol does not establish clinical utility, causal identification,
calibrated uncertainty, cell-response distributions, genetic interactions,
cross-context biological prediction, Bayesian/mechanistic coupling or GPU speed.
