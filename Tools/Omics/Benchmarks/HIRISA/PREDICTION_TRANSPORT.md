# HIRISA aggregate prediction transport supplement

2026-09-10, after the synthetic batch product gate and before any HIRISA
response-model fitting. This supplements, and does not replace, the original
`PREDICTION_EXECUTION.md` and its 79-fold manifest. DE outcomes have already been
inspected. No model settings or fold selections depend on those outcomes.

The new `singlecell-perturbation-batch` API reuses one native pseudobulk bundle.
It copies the original H5AD, mapping, report and publisher receipt into private
staging and reconstructs the entire report from original counts once per
publication or verification. APFS cloning can share physical extents; the files
remain independent snapshots. The old publisher receipt is retained unchanged.
The current executable executes reconstruction even when the old publisher was
a different build; the batch receipt binds the current implementation as well
as the original publisher receipt and exact source/report identities.

Each fold explicitly selects the frozen training library group indices and only
the held-out control group. Its population becomes the explicit aggregate cell
group. Original library IDs, donors, batches, integer counts, gene order and
global source-cell indices remain unchanged. Donor, sample and source-cell
overlap are rejected before fitting. The fitting and prediction primitives
receive isolated aggregates, with no held-out treated rows or target URLs.
The scorer's target accession stays in a separate scoring manifest.

`VivoPerturbation.swift` now delegates its existing report-based methods to
aggregate overloads. Normalization, context selection, scaling, alpha, solver,
clipping and all four baselines retain the frozen numerical implementation.
`VivoPerturbationIO.swift` is unchanged. The transport receipt records both the
old frozen owner hashes and the new source hashes. Model source fingerprints
identify the isolated training aggregate; the outer receipt supplies raw-H5AD
provenance. Altering held-out treated counts therefore cannot alter fitted
model bytes or predictions.

Seven focused Swift tests and ten synthetic product checks passed before this
supplement. The latter include independent count/membership and NumPy model
oracles, cross-build raw-source reconstruction, byte-exact repetition, and
rejection of rehashed models, rehashed source counts, overwritten outputs and
fabricated outputs in failed folds. These qualify software boundaries only.

The real experiment must retain all 79 fold statuses and full model/query
outputs. Freeze and hash every output before the scorer opens treated counts.
Then independently reconstruct every fit from training/control counts, and
score both predeclared gene families without parameter changes. Failed folds
make the corresponding primary contrast summary incomplete; do not average
only successful folds. Native reconstruction and exact replay remain separate
from empirical held-out prediction performance and biological validity.

This transport avoids 79 repeated whole-source scans. It still retains the
source report and cell metadata in memory and is not a claim of fully
out-of-core metadata or GPU acceleration.
