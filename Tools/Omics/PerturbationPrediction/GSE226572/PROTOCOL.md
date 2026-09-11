# Complete GSE226572 PBMC response-transfer experiment

Declared 2026-09-11 after study design and the complete primary SOFT metadata
were read, before downloading or inspecting any count matrix in this cohort.
Prior Kang/HIRISA prediction failures and annotation-transfer diagnostics have
been inspected. This is a new prediction evaluation, not blinded model discovery.
It leaves the separate GSE181897 experimental-code admission unresolved.

## Endpoint and primary experimental identities

Use all 24 original raw 10x H5 matrices listed in the GSE226572 family SOFT.
The source identifies three healthy male donors, 19–22 years old; retain literal
sample-title donor identifiers D34, D38 and D39. Time fields explicitly identify
0, 1, 2, 4, 8, 12, 24 and 36 hours, with donor-specific available times. Preserve
the six separate zero-hour libraries, two per donor, and all 18 stimulated
libraries. Do not invent missing donor/time combinations or count libraries,
time points, genes or cells as independent biological donors.

The primary extraction protocol explicitly states recombinant IFN-beta 1a at
1000 U/mL final concentration, staggered addition and collection after 36 hours
of culture for every sample, including unstimulated controls. Thus zero-hour
libraries are culture-matched unstimulated controls, not fresh cells. Source
sample titles and time fields must agree; retain every source accession and URL.
Do not infer interventions from expression. Genome annotation is GRCh38 release98.

The endpoint is the aggregate RNA profile of all cells passing the source's
stated initial QC: at least750 detected genes and mitochondrial RNA fraction
at most20%. Compute these exactly from all measured RNA counts, with positive
counts defining detection and exact source symbols starting MT- defining the
mitochondrial set. Retain every barcode's counts, detection, mitochondrial total
and inclusion mask. No upper gene/count cutoff, response-dependent filter,
cluster removal, doublet classifier, or learned identity is introduced. This
reproduces only the stated initial count-QC criteria, not the authors' subsequent
cluster exclusions or115503-cell curated object. Do not call the admitted
barcodes authoritative PBMC identities or qualified singlets.

This whole-population endpoint is intentionally distinct from prior B-cell
studies. It includes composition, viability, preparation and library effects;
it does not isolate cell-intrinsic regulation. The processed9.7GB Seurat object
is not needed to define this endpoint. No unseen source label is invented.
Aggregate the two zero-hour libraries of each donor by summing their counts;
report both source libraries separately in QC and provenance. Normalize the
pooled aggregate by its full RNA total. Keep every admitted stimulated library.

## Frozen model, panel and evaluation

Train the existing native no-change, mean, median and context-ridge baselines on
all eight Kang donors' complete admitted PBMC populations and both ctrl/stim
conditions. Use the source-qualified full Kang count file, with all source cell
types included and each source library's full measured-gene denominator retained.
The Kang training source has been inspected previously. The query cohort is
held out in its entirety: no model selection, query response fitting, scaling,
calibration, feature selection by query expression or threshold tuning.

Use exact unique gene-symbol matches between the full Kang feature set and every
query source RNA feature set. Retain original Ensembl identifiers and symbols,
all ambiguous/unmatched features and reasons. Symbols do not establish a common
annotation release. Do not zero-pad missing features; all nonpanel RNA features
remain in normalization denominators. Require all24 query files to agree on
feature identities, or explicitly resolve source metadata before fitting.

Keep existing model defaults: log1p(CPM), context-ridge alpha1, at least10 training
counts and at least2 expressing training donor pairs, training-only scaling,
mean-response intercept and nonnegative clipping of predicted treated values.
The only query input is each donor's observed pooled zero-hour profile. Freeze
all three native query predictions before comparing with any stimulated profile.
The model has no duration covariate: its transferred six-hour Kang response is
evaluated unchanged against every available query time. This explicitly tests
how far that fixed response transfers, not a learned temporal response model.
Dose, disease status, culture duration, assay and cell composition also differ.

Primary criterion: the native mean-response method must reduce all-panel RMSE
by at least5% versus no change when averaging first over all available stimulated
times within each donor, then equally across all three donors. Report every
individual time point and donor, per-time means/support, RMSE, MAE, response
correlation, negative clipping and worse-than-no-change counts for every method.
The ridge method's secondary comparison must beat both no change and the mean.
Do not select a favorable time, donor, gene subset or method after scoring.

Report nominal95% future-donor mean-response intervals if supported by the
existing native owner, using training donors only. Retain pointwise coverage,
miss directions, widths and availability for every donor/time; do not call them
calibrated, simultaneous or a temporal model. Eight training and three query
donors do not provide broad population or clinical qualification. No endpoint
here establishes protein function, tissue behavior, disease or treatment benefit.

## Complete execution and storage

Validate every raw sparse record, shape, index ordering, count domain, feature
and barcode identity using bounded chunks; never densify cells by genes. Retain
all source files with byte counts and SHA256. External preparation into H5AD is
separate from native import, aggregation, model fitting and prediction. Compare
every selected cell/feature/count with the original source and every native
aggregate and group membership with independent sparse arithmetic.

Freeze this protocol and all source metadata before counts; freeze count inputs,
feature mapping and model plans before fitting; freeze all predictions before
scoring. Reconstruct native bundles, independently verify all model outputs and
repeat scoring for identical results. Preserve failures and complete-cohort
counts. Use the laptop because the Mac mini lacks capacity; do not duplicate
source files there. Archive new evidence and remove working copies only after
hash checks, successful restoration and open-handle checks.
