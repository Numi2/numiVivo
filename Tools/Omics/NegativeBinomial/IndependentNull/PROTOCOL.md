# Independent-donor randomized null benchmark

Frozen 2026-09-10 before count aggregation, model fitting or any DE outcomes.
Base implementation: 73aea31c00a7adfc9743b4373efa3f2c17e73b08.
Only source metadata and matrix storage dimensions have been inspected.

Use the complete public CELLxGENE Human Immune Health Atlas B and Plasma cell
release, dataset ed78c6b7-cb40-4ebc-b55f-40bc99ba1b5d, version
b6986a7f-981e-4e04-93ed-57f444749b8b. Source SHA256:
2f019ba5ae48fdf61e02f510aac66cf6035b44aa4c582d52efce9ebccad90ee3.
Keep all 160,632 observations and all 32,357 raw-axis features. Read integer
counts from raw/X, validating every stored value; never round the float storage
or substitute normalized X. Preserve all original observation and feature IDs.
The author AIFI_L1 label is B cell for the entire released matrix, including
its five AIFI_L2 populations. Use that complete broad population. Retain subtype
composition descriptively; do not select a subtype after observing DE outcomes.

The release is a curated, processed representation of observed UMI counts,
marked is_primary_data=false because these cells also occur elsewhere in the
CELLxGENE corpus. Do not combine those other releases as independent observations.
The 108 donor_id values identify different people. Require one original sample
and one reported batch per donor. Preserve the original demographic/CMV and
sample provenance, without treating age or infection status as a sham treatment.

Create nine disjoint cohorts of twelve original donors. Sort donors first by
literal original batch_id, then SHA256 of UTF-8
`numivivo-independent-null-v1|cohort|donorID`, then donorID; consecutive groups
of twelve form cohorts 01 through 09. No donor or cell appears in two cohorts.
Within each cohort and source batch, sort donors by SHA256 of
`numivivo-independent-null-v1|arm|donorID`, then donorID. Assign alternating
positions to shamA and shamB. For even-sized batches start at shamA. Sort odd-
sized batches by SHA256 of `numivivo-independent-null-v1|odd|cohortID|batchID`
and alternate their starting arm shamA/shamB. Since each cohort has twelve
donors, the number of odd batches is even and each cohort has six donors per arm.
Assignment uses no count, feature, effect, test probability or cell-type fraction.

Use the existing native streamed H5AD publication/replay owner, preserving the
complete source bytes. Independently verify every cell QC value and all raw
pseudobulk counts against chunked SciPy aggregation. No dense cell × gene matrix.
The current 1 GiB source-byte admission rejects this 1,199,390,735-byte source;
repair that structural limit using the existing fixed-buffer snapshot reader,
without reducing the source or changing count/design/fit budgets to fit outcomes.
Retain the original rejection and measure native publication/replay memory.

For each cohort run native NB Wald, LRT and modern adjusted QL on the existing
independentReplicates design with original batch fixed effects, native median-
ratio offsets, gammaParametric dispersion trend, minimum twenty trend genes,
minimum prior variance 0.25 and outlier threshold two standard deviations.
Use at least three donors per condition, ten cells per pseudobulk, total gene
count ten and expression in three pseudobulks. No influence exclusions, count
replacement, active-donor policy, effect prior, trend fallback or seed selection.
Use the original includedDonorIDs contract; preserve every withheld feature and
failed cohort. Do not replace a failed cohort with a smaller or easier one.

Use pinned edgeR 4.10.5 modern robust QL, limma 3.68.5 voom/robust empirical
Bayes, and the existing pinned DESeq2 Wald reference. All consume exact original
aggregate counts, the same numeric batch design and native offsets. Preserve
reference-specific gene families, warnings, failures and DESeq2 filtering; also
report explicit unfiltered BH comparisons. Do not tune solver/trend choices
against the sham outcomes. Record any numerical sensitivity as a separately
identified follow-up rather than silently replacing the declared result.

Report attempted/eligible/tested/withheld genes, fit failures, dispersion
boundaries, raw P fractions below 0.01/0.05, BH calls at 0.01/0.05/0.10, and any
BH ≤ 0.05 event for every cohort/method. Report full method families and a
separate jointly tested family. Keep all probabilities/effects and count/design
identities. Check native score/likelihood/F-tail/BH arithmetic independently;
these numerical checks must not be described as calibration tests.

These condition labels are randomized and cause no biological intervention;
discoveries are false under the assignment null. Cohorts use disjoint donors,
improving on overlapping within-library cell splits, but remain nine experiments
within one selected study with shared collection/processing practices and broad
B-cell mixtures. They are not nine independent studies or a precise universal
FDR estimate. Report the nine-cohort event fraction descriptively without an
IID binomial confidence claim. No alternative-model power or effect-interval
coverage is tested. This evidence alone cannot select a production default.

Sources: [author data description](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/scrna/),
[author cohorts](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/cohorts/),
[CELLxGENE collection](https://cellxgene.cziscience.com/collections/77f9d7e9-5675-49c3-abed-ce02f39eef1b),
and [Gong et al.](https://doi.org/10.1038/s41586-025-09686-5).
The randomized null design is ours, not a reproduction of the paper's aging analysis.
