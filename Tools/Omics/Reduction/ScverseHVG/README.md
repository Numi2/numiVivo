# Full-cohort highly-variable-feature comparison

Executed 2026-09-12 against the same Scanpy 1.12.4 environment used for the
[matched PCA comparison](../ScversePCA/README.md). This independently selects features
from the original sparse counts; it does not feed native selections into Scanpy.

Each full source matrix is converted to float64 sparse storage, normalized to 10,000
counts per cell across all source features, then log1p transformed. Scanpy runs
`highly_variable_genes(flavor="seurat", n_top_genes=2000, n_bins=20)` without a batch
key. Every source feature identity and its index are checked against the native
model before comparing the selected index sets. No dense cells × genes matrix is
constructed; this reference keeps the sparse matrix resident.

| Cohort | Cells | Source features | Native/Scanpy selected | Shared |
| --- | ---: | ---: | ---: | ---: |
| Kang | 24,673 | 15,706 | 2,000 / 2,000 | 2,000 |
| Hagai | 13,863 | 22,048 | 2,000 / 2,000 | 2,000 |

Both pass the pre-run exact-membership gate, with empty native-only and Scanpy-only
sets. Maximum finite log-mean differences are 1.78e-14 and 1.11e-14; log-dispersion
differences are 8.30e-14 and 1.07e-13. Normalized dispersion differs by up to 9.41e-7
and 8.97e-7, respectively, without changing selection. Hagai has 6,599 nonfinite
log-dispersion and normalized-dispersion entries in each implementation; these are
retained in the full feature table, not converted to successful finite comparisons.

Together with the separate matched-feature PCA run, this establishes that these
two full real cohorts produce the same selected feature membership and numerically
matching PCA results. It does not establish all Scanpy flavors, tied-cutoff behavior,
batch-aware feature selection, million-cell execution, or biological generalization.

The one-shot HVG stage times were 0.531 s Kang and 0.123 s Hagai. They exclude
input, normalization, imports and PCA and are not a native/scverse speed comparison.
`run.py`, `protocol.json`, results and log retain the actual execution. Compressed
CSV files preserve every Scanpy feature statistic. `manifest.json` hashes the raw
and compressed tables and execution files; original native model hashes are in
`results.json`, bound to the earlier full PCA evidence. Runtime files remain at
`/Users/n/numivivo-hvg-scverse-20260912`. The driver uses the same single-thread
numeric environment as the PCA comparison and requires the original H5AD/model
paths. Reproduction should use a fresh output root.
