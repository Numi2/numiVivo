# Frozen HIRISA inference execution

2026-09-10, before any HIRISA model fitting, at source commit
981e3f456fd63a143c218e190103d34a551a63d2. The original PROTOCOL.md remains
unchanged. Full source and backed AnnData checks have completed. Complete-source
native ingestion is still queued behind the active transfer.

Prepare each of the sixteen predeclared donor-matched contrasts from the exact
per-library source-audit gene counts, which already agree with the complete
AnnData source. Include every cell UUID and every gene in each selected library;
record the original whole-source observation ranges and per-library reference
hashes. Group these observations by their experimentally supplied enrichment,
keeping each accession as a separate pseudobulk row and the original donor as
the biological replicate. This is an explicit Python-prepared count input to
the native inference owner. It does not claim native raw-count aggregation.
After full native ingestion, require exact agreement of all library counts
before connecting the native raw-input and inference evidence chains.

Run the existing shared VivoPseudobulkDifferentialExpression.evaluate owner via
the previously compiled, source-hash-verified IndependentNull harness. The
harness's historical executable name does not assign null labels here: inputs
retain the original experimental treatments. Do not reimplement native fitting
in Python or change its numerical settings based on these outcomes.

Use all 48 protocol requests: sixteen Wald, sixteen LRT and sixteen adjusted QL.
Run all Wald cases first, then LRT, then QL, with stable population/treatment
ordering from the frozen design. Keep every terminal error and failed fit.
The control is none; the treatment is each original IFN label. Use pairedDonors,
medianRatio, adjustForBatch=false (batch is matched within donor and aliased with
donor terms), minimum three donor pairs, ten cells per pseudobulk, ten total
feature counts, expression in three pseudobulks, and existing minimum ten
all-positive normalization features. Use gammaParametric, minimum twenty trend
genes, minimum prior variance 0.25 and dispersion outlier threshold two SD.
No influence threshold, active-donor policy, effect prior or fallback is added.
The four-pair Bcell IFNg case keeps exactly the predeclared paired donors.

Prepare the exact numeric donor design and all-positive median-ratio offsets
independently from these aggregates. Verify them against native results. Use
the existing pinned edgeR/limma/DESeq2 reference harness without method changes;
its numeric design input supports these donor terms. Preserve reference warnings,
errors, default DESeq2 filtering and unfiltered BH results. Compare full families
and a separately identified jointly tested family. Do not describe correlated
methods as independent truth, or this study as FDR/power/coverage calibration.

Run count/design/likelihood/score/information/tail/BH checks against the native
outputs. Treat numerical consistency, software interoperability, experimental
response concordance and predictive validity as separate evidence. Prediction
settings remain unfitted and require their separate pre-fit specification.
