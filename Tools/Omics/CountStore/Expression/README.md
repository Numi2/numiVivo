# File-backed counts to native differential expression

The native count bundle now feeds the existing negative-binomial differential-expression model without reconstructing resident cell-membership lists. This closes the count-to-analysis interface for the complete Parse PBS/IFN-beta cohort: **725,031 cells, 40,352 features, 12 donor pairs and 24 aggregate observations**.

This estimates differences between **observed** pooled PBMC responses. It does not predict an unseen response or distinguish cell-intrinsic regulation from changes in cell composition. The [biological prediction assessment](../../../../Documentation/BiologicalPrediction.md) retains the independent-validation failures and unresolved Parse prediction inputs.

## Interface and ownership

```sh
numivivo-omics singlecell-file-expression counts-bundle --plan contrast.json --output expression-bundle
numivivo-omics singlecell-file-expression-verify expression-bundle
```

The schema-1 plan binds the exact source receipt and source implementation separately from the analysis implementation, plus the existing contrast. See the retained `file-plan.json` for the executed plan. The analysis owns a private immutable source snapshot. Opening verifies receipt bytes; execution scans every QC/membership row and checks group sizes and aggregate library totals before calling the existing statistical core. Verification recomputes the statistical report. Reconstructing raw feature counts still requires the count owner's original-stream verification.

`VivoOmicsDesignObservation` represents either original resident indices or a receipt-bound file membership reference. File observations expose the actual `sourceCellCount`; `sourceCellIndices` is absent, never a fabricated empty list. Legacy observation JSON remains byte-compatible. Minimum-cell filtering, donor pairing, design construction and Wald/LRT/QL arithmetic share the same owners. Swift callers now receive `[VivoOmicsDesignObservation]` from `design.observations`; use `sourceCellCount` for group size and handle optional `sourceCellIndices` when accessing membership. This is a source API change even though resident JSON is preserved.

## Actual qualification, 2026-09-11

The 107-input scoped native build and **35 tests across six suites** pass, including exact legacy encoding, Wald/LRT/QL agreement, empty groups, distinct source/analysis identities, source mutation isolation and rejection of a rehashed incorrect statistical report.

Both full-cohort native fits and statistical reconstruction pass. Every non-design JSON value and every design value except the explicit membership representation matches the prior native owner exactly. The baseline helper calls the actual old native library; it is a native-core runner, not the old product's entire command workflow.

The frozen paired NB design uses median-ratio size factors, parametric dispersion trend with no fallback, 12 donor pairs, 13 columns and 11 residual degrees of freedom. Of 40,352 source features, **33,899 are tested, 5,946 are rejected for rank-deficient support, and 507 fail the low-expression filter**. All statuses and feature rows are retained.

Six real R workflows completed with edgeR 4.10.5, limma 3.68.5 and DESeq2 1.52.0 under R 4.6.1. They use the same complete aggregates, count filter and paired design. Shared native offsets and each package's own normalization are reported separately. Independent design, libraries and size-factor checks pass. All six runs completed without recorded warnings or messages.

| Workflow | Effect RMSE | Effect Spearman | P-value Spearman | BH 0.05 set Jaccard |
| --- | ---: | ---: | ---: | ---: |
| Shared native offsets: edgeR-QL | 0.076546 | 0.994101 | 0.987019 | 0.878102 |
| Shared native offsets: limma-voom | 0.359391 | 0.859523 | 0.884023 | 0.625522 |
| Shared native offsets: DESeq2 | 0.134987 | 0.989723 | 0.995313 | 0.916705 |
| Package normalization: edgeR-QL | 0.101727 | 0.993925 | 0.880909 | 0.762155 |
| Package normalization: limma-voom | 0.336878 | 0.859670 | 0.818986 | 0.575170 |
| Package normalization: DESeq2 | 0.134954 | 0.989734 | 0.995314 | 0.916563 |

Comparisons use 33,899 finite paired effects, within 39,845 count-filter-eligible features. There are 15,167 native BH-0.05 calls in each comparison; the complete reference call counts, discrepancies and exclusions are in [reference-comparison.json](evidence/2026-09-11/reference-comparison.json). QL F, moderated-t and Wald tests need not return identical p-values. This is descriptive concordance, not measured false-discovery calibration or biological ground truth.

The file-backed command reached **970.4 MiB native peak RSS in 376.2 seconds**, versus **984.4 MiB in 373.3 seconds** for the prior native-core runner. The complete feature diagnostics and JSON report still impose substantial resident memory. These results do not support an analysis-memory or speed advantage. The subsequent [incremental report writer](Memory/README.md) preserves every report byte and lowers native peak RSS to 202.3 MiB; the timings and memory above remain the original implementation comparison. Native measurements include loading, inference and writing; the baseline's Python gzip transport is outside native metrics. R method timings have a narrower scope and are not an end-to-end speed comparison.

## Evidence and reproduction

[The evidence manifest](evidence/2026-09-11/manifest.json) binds the complete native analysis bundle and its count snapshot, both native executables, frozen native source, plans, complete R inputs/outputs, feature comparison table and logs. The source count validation is a separate [cell-axis/count qualification](../CellAxis/README.md). No whole 227 GB source-file hash or new prediction outcome is claimed. Parse data and derivatives are **CC BY-NC 4.0**.

The prior native report is retained losslessly as the new report plus its exact original observation block. Preparation directly compared all 113,012,281 reconstructed native JSON bytes to the original output (SHA-256 `83be008ef62c9b878f14afb3f52cfdc540011e8a9de1130a90f5bfd50a8b0883`). The original Python gzip wrapper is not reproduced. This saves duplicate storage without replacing measured output with a new model run.

`retain.py verify --archive evidence/2026-09-11/results.tar.gz` checks every archived object. `retain.py restore --archive ... --out fresh-directory` extracts verified bytes; optional `--study original-study` allows hash-checked separate-inode APFS copies to conserve disk. Then use the restored executable to run `singlecell-file-expression-verify fresh-directory/native-file`, and `report_delta.py verify --study fresh-directory` to reconstruct the baseline digest. Restoring files and rerunning statistical inference are separate checks. Both now pass in the [retained restoration evidence](evidence/2026-09-11-restoration/manifest.json): all 186 archived members verify (67 source APFS copies and 119 extractions), the exact baseline JSON reconstructs, and the restored native executable independently recomputes the complete statistical report with exit 0.

The [incremental writer](Memory/README.md) now removes whole-report JSON materialization. The next memory targets are resident per-feature fit/diagnostic state and broader observation/feature scaling. The next prediction target remains a predeclared independent experiment with resolved feature identities, dose and outcome definition; this DE comparison does not waive those requirements.
