# Replogle 2020 UPR: independent fixed-predictor validation

Declared 2026-09-11 at owner base `50ab30168c5dddd5ee338249365a2bb361ab768e`, after inspecting the paper and original guide/barcode metadata but before downloading the expression matrix, fitting or scoring this study. Adamson's existing frozen protocol and unresolved identities remain unchanged.

## Source and inclusion

Use the complete original deposited experiments 1–5 in GSE146194 (GSM4367979–GSM4367983), Replogle et al., Nature Biotechnology 2020, DOI 10.1038/s41587-020-0470-y, PMID 32231336. GEO supplies the shared exp1-5 Matrix Market count matrix, features, barcodes and cell identities through GSM4367979. Retain all 40,997 source barcodes, source features, counts and assignment records. GEO identifies UMI-collapsed Cell Ranger 3.0.0 counts aligned to GRCh38-1.2.0. Preserve full barcode suffixes and explicit gemgroups; never join truncated barcodes or infer platforms from expression.

Audit all five gemgroups separately. Verify the gemgroup-to-platform mapping against author metadata before naming them as platforms in predictive results. They are technical/capture conditions in K562, not independent donors. No differential-expression significance or biological replication is inferred from these groups.

The primary paper's Pilot UPR methods explicitly identify sgNegCtrl2 and sgNegCtrl3 as the two non-targeting controls. Keep both counts separately in the audit; pool their raw counts within each platform only when forming that platform's prediction control. Use the deposited good_coverage and number_of_cells fields to retain uniquely assigned, confident cells: good_coverage=True and number_of_cells=1. Retain every excluded row and reason; require unique full-barcode joins. No expression-effect, gene-variance, cell-library-size or outcome-score threshold selects the cohort, except that any zero-library condition cannot enter a logarithmic model and must be reported.

Retain all 30 non-control guide groups. Experimental guide identities come from original author assignments. Gene-level descriptor mapping requires explicit, unambiguous identity evidence; a matching expression symbol alone is insufficient. Do not guess aliases or force unsupported guides into a gene. Keep unresolved mappings and generic-baseline coverage explicit.

## Frozen predictor and endpoints

Use the existing native VivoTargetKernel unchanged: full-source-gene log1p CPM denominators, direct BP/MF/CC term Jaccard, lambda 1, unpenalized intercept, clipping at zero without CPM reclosure. Use the same exact-identity MyGene/GO capture rules as the frozen Adamson protocol, excluding NOT and GO roots with no ancestry expansion. Record source build and annotation citation overlap with this and Adamson studies; current knowledge is not temporally independent.

For every admitted target and each of five technical platforms, train within that platform on controls and all other unambiguous single targets. Withhold every condition for the query target. Fitters receive only selected training counts and descriptors; query inputs contain no held-out expression. Freeze exact mappings, all folds and prediction payloads before the separate scorer reads held-out response vectors. No hyperparameter selection or outcome-driven tuning.

Primary endpoint within each platform: equal-target mean all-gene log-response RMSE on the same supported targets, against all-training-target mean, supported-training-target mean and the fixed one-position supported-response shuffle. A useful GO-specific gain must beat all three. Report all five platform outcomes, every fold, no-change regressions and coverage; no favorable-platform selection or pooled platform pseudo-replication. Secondary training-selected top-1000 mean-CPM results, MAE, correlation and sign metrics cannot rescue failed primary results.

This independently collected study tests the fixed algorithm; it does not transfer Norman-trained coefficients. The study reuses UPR targets selected from Adamson and shares investigators and experimental systems. It does not establish unseen biological contexts, prospective target selection, uncertainty coverage, new tissues or clinical prediction.

## Verification

Before fitting, compare every converted H5AD count and axis with the original Matrix Market/GEO files, then execute native streamed pseudobulk aggregation and reconstruction. Compare every group sum, membership and cell QC with an independent sparse reference. Never allocate a dense cells-by-genes array. Verify target-level training exclusion by mutating withheld counts. Bind exact source, scripts, executable and OS identities; compare all native prediction vectors with the existing independent centered-kernel reference. Preserve failed and incomplete checks. An ingestion or metadata pass is not prediction qualification.

## Assay clarification before count conversion or fitting

The original feature TSV contains 33,694 Gene Expression and 64 CRISPR Guide Capture entries. The first conversion admission correctly rejected the assumption that all features were RNA. Retain all source entries and exact feature-type metadata, then partition both assay feature axes with all cells retained. Only original Gene Expression features enter RNA denominators and prediction panels. This source-declared assay distinction is not outcome-based feature selection. Full mixed-source and RNA-specific totals/nonzeros are reported separately.
