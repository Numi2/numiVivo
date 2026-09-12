# Explicit identities with numeric dataframe indices

Native annotation, count-store import and pseudobulk can now operate on numeric-index H5AD files when the import plan explicitly selects string identity columns. The count importer validates axis lengths independently, then reads only the identity columns required by the plan. Previously, it decoded both original indices even when every relevant identity had an explicit override.

Set all three of `barcodeColumn`, `featureIDColumn` and `featureNameColumn` to the appropriate source columns. `sampleColumn` remains required and sample/design definitions remain explicit. If a mapping is absent, the historical index fallback remains in force: numeric indices are not silently stringified or promoted to biological identifiers. Missing selected identities, duplicate cell identities and dimension disagreements still reject publication. For `raw/X`, feature columns must exist in `raw/var`; annotation does not invent raw feature mappings.

Native annotation accepts numeric-array and nonmissing nullable integer/boolean indices while retaining their values. Existing rejection of missing nullable dataframe indices during annotation remains. The count importer can validate unused nullable index shapes independently of the explicitly selected identifiers; these are distinct operations.

## Execution evidence

- 24 complete annotation → count-store → pseudobulk routes cover `X` and `raw/X` across twelve legacy/modern numeric-index fixtures. Source datasets and dataframe semantics survive annotation, every nonzero count matches an independent sparse AnnData reference, explicit barcode/feature identities match, all aggregate gene sums agree, and count/aggregate replay passes.
- Seven rejection checks cover omission of each identity override, short barcode columns, duplicated or missing selected barcodes, and malformed index rank. No rejected output is published.
- The new binary reconstructs the complete retained real Kang count store: 24,673 cells, 15,706 genes and 14,184,532 count records. Native reconstruction recomputes and checks the original metadata, QC and count hashes. It also reconstructs the original 124 pseudobulk groups and 866,103 aggregate nonzeros against the retained report. This is actual source execution through the updated reader, not report parsing alone.

The conformance fixtures are synthetic and their numeric indices are not inferred gene identifiers. The Kang check is a real-data software regression, not new biological validation or a new prediction score. This scoped build does not substitute for a full product build, H5MU-wide identity support or unrestricted numeric identity coercion.

[Retained evidence](evidence/2026-09-12-explicit-identities/manifest.json) includes the earlier annotation/count failures, all source inputs, executable, checks, seven rejection diagnostics and Kang reconstruction logs. The initial compile failure and a corrected harness expectation about CLI output are also retained. Large original Kang bundles remain separately retained and hash-bound in `kang-regression.json`.

Run `check_explicit_identities.py --binary H5AD_CHECK --fixtures NUMERIC_INDEX_CHECKS --out NEW_DIRECTORY`. The fixture directory is the `checks` directory from the preceding numeric-index evidence archive. Native annotation and count execution require no Python; AnnData/NumPy are used only by the independent checker.
