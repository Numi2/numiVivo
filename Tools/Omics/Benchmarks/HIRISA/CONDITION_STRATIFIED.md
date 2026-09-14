# Condition-stratified integration candidate

The opt-in condition-stratified donor correction was executed on the complete
HIRISA PCA bundle on 2026-09-14. It fits the categorical donor effect separately
inside each declared treatment condition, using the existing fixed ridge and
the existing soft-clustering objective. The historical global correction remains
the default; adaptive penalties and the joint donor-plus-batch route are rejected
for this scope by plan validation.

The native publish and replay both passed on the M4 Pro. The run processed all
**1,612,594 cells**, 20 PCs, 100 clusters and five donor levels. The candidate
stopped on relative objective tolerance after six correction iterations; its
maximum ridge residual was `5.809779771020291e-14`. The output contains the
complete scores, memberships, assignment scores, metadata, report and receipt.
The production verifier reconstructed the parent PCA and the integration output
before accepting the receipt. A pinned native HDF5 library was supplied through
`NUMIVIVO_HDF5_LIBRARY`.

## Preservation result

The frozen donor-held-out diagnostic keeps all 114 folds, 471 observed
stratum/label comparisons, 330 metadata-defined rare comparisons and the same
146 supported, control-sensitive comparisons used by the existing native result.
Author `celltype.l2` values are comparison annotations, not ground truth. The
diagnostic is transductive because the full-cohort PCA and integration precede
the donor-held-out decoder; it is not a prospective identity or outcome model.

| Comparison | Supported / sensitive | Margin failures | Rare failures | Mean recall loss | Complete gate |
| --- | ---: | ---: | ---: | ---: | --- |
| Original PCA → native condition-stratified candidate | 146 | 40 | 10 / 29 | -0.00749324 | **fail** |
| Existing native global correction → candidate delta | 146 | 11 | 3 / 29 | -0.00377681 | **fail** |
| Existing native global correction | 146 | 37 | 10 / 29 | — | **fail** |

Negative mean loss means the candidate's equal-donor mean recall was higher. The
aggregate improvement does not satisfy the complete gate: the candidate has 40
absolute margin failures, three more than the existing native correction, and a
maximum individual-donor loss of `0.38089786`. The candidate is therefore
**rejected for promotion**. Its result is useful as a measured local-alignment
experiment, not evidence that condition labels, cell types or biology were
preserved.

## Numerical and artifact checks

The independent per-cell augmented-SVD oracle passed all 114 folds for both the
absolute and native-relative result files. Maximum coefficient error was
`6.869504964868156e-16`; maximum confusion error was zero. All candidate output
hashes, the plan, executable, HDF5 dependency and failed admission attempts are
bound in the [compact evidence manifest](evidence/2026-09-14-condition-stratified/manifest.json).
The two failed attempts are retained: the unmodified old PCA receipt was
rejected for executable identity, and the first regenerated-input invocation
was rejected until native HDF5 was supplied.

The full 516,030,080-byte score matrix, 2,580,150,400-byte memberships matrix,
assignment scores and source H5AD remain external on the Mac mini under the
hashes recorded in `evidence/2026-09-14-condition-stratified/execution.json`.
The compressed [absolute result](evidence/2026-09-14-condition-stratified/condition-f1-vs-baseline.json.gz)
and [native-relative result](evidence/2026-09-14-condition-stratified/condition-f1-vs-native.json.gz)
retain every fold, confusion matrix and measured comparison.

## Reproduce

Use the published `f1b882f1b4d545ab8a40e149a514e4a1b2be8a99` binary on the Mac
mini, the frozen HIRISA PCA source and a new output directory:

```sh
export NUMIVIVO_HDF5_LIBRARY=/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib
export DYLD_LIBRARY_PATH=/Users/n/numivivo-multiassay-hdf5-20260909:/opt/homebrew/opt/hdf5/lib
numivivo singlecell-pca-integrate \
  /Users/n/numivivo-hirisa-20260910/integration-condition-20260914/input-candidate \
  --plan /Users/n/numivivo-hirisa-20260910/integration-condition-20260914/plan.json \
  --output /new/condition-stratified/native
numivivo singlecell-pca-integrate-verify /new/condition-stratified/native
```

Run `check_annotation_retention.py` with the frozen `annotation-retention`
`rows.npz`, `ledger.json` and `freeze.json`, once against `baseline.json` and
once against `native.json`; then run `check_annotation_retention_oracle.py`
for each result. Never overwrite the completed roots. These checks establish
artifact and numerical consistency only. Independent biological endpoints,
local composition/alignment, unseen contexts and phenotype outcomes remain open.
