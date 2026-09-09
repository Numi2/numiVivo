# Full Norman 2019 perturbation benchmark input

This extends the native streaming input path to the complete published
**filtered** scPerturb Norman 2019 release. It prepares the next perturbation
prediction benchmark; it does not establish a fitted unseen-perturbation model.

The source is [scPerturb RNA release 1.4, Zenodo 13350497](https://zenodo.org/records/13350497),
file `NormanWeissman2019_filtered.h5ad`, from
[Norman et al., Science 2019](https://doi.org/10.1126/science.aax4438),
GEO GSE133344. The complete filtered release contains 111,445 K562 CRISPRa cells,
33,694 features, 361,582,621 stored count entries and 1,635,387,239 UMIs.
It has 105 single-target conditions, 131 paired-target conditions and control.
This is the published filtered matrix, not the original experiment's raw release.

## Source and biological identities

`prepare.py` pins the exact byte size, published MD5 and independently computed
SHA-256 before reading. It checks every CSC coordinate and count: sorted indices,
no duplicates, no zero/negative/nonfinite/fractional values, and unique cell and
feature identities. It keeps all cells and all features. The source stores
counts as float32, but every value is an exactly representable nonnegative integer.

The source's `ensemble_id` column (that spelling is in the file) supplies the
Ensembl feature IDs; its gene-symbol index supplies feature names. Mitochondrial
QC uses source symbols beginning `MT-` with their matching Ensembl IDs. The
mapping retains the original AnnData file byte for byte.

All observed controls carry explicit negative-control guide pairs. The upstream
scPerturb processing code can also classify `no_reads_found` as control; the
preparation rejects that possibility for this pinned release instead of silently
mixing missing-guide cells into controls.

The eight `gemgroup` values are recorded in the audit. They do **not** establish
eight biological replicates or donors. Condition aggregation uses one explicitly
unresolved pooled biological-replication identity, no donor ID, and a pooled
batch label. No DE contrasts or inferential replicate tests are requested.
The original source retains per-cell guide and GEM-group metadata for later
batch-aware analyses.

The descriptive activation audit directly matches 102 of 105 single-target
names to source symbols. `C19orf26`, `C3orf72` and `KIAA1804` need explicit alias
resolution; no gene ID is guessed. Of the 102 matched targets, 100 have a positive
pooled target-gene log1p-CPM change. `CELF2` and `CKS1B` do not. These observations
are retained, not treated as replicate-level DE or proof of causal activation.

## Evaluation splits

The prepared `splits.json` contains two distinct categories:

* **Unseen combinations:** train on the 105 singles, with control available;
  hold out all 131 pairs together. Every pair has both constituent single targets
  represented in training. Pair outcomes must remain sealed until scoring.
* **Unseen targets:** 105 separate target folds. Each removes every single and
  paired condition containing the held-out gene from training. Merely withholding
  one pair is not an unseen-target test. A model needs independently supported
  target descriptors to make target-specific predictions in these folds;
  a one-hot unknown-gene fallback is not such evidence.

These are condition-mean split definitions, not validated predictions or
single-cell response distributions. They do not demonstrate donor generalization,
and they do not connect expression changes to Bayesian mechanistic parameters.

## Native resource change

The shared streaming source allowance is now one billion stored entries. Source
bytes remain capped at 1 GiB, aggregate nonzeros at five million, and input reads
at 65,536 entries. A sparse major segment, metadata, per-cell QC, grouped sums and
the serialized report remain resident. The resident matrix import and general
H5AD axis-projection limits are unchanged.

The old product rejected this exact source at its 100-million-entry guard; the
failure is retained as `prior-product-rejection.log`. The full source aggregates
to 3,577,278 nonzero condition-gene values, within the existing aggregate bound.
The new oversize regression uses unallocated chunked arrays declaring
1,000,000,001 entries, so it checks rejection before a large scan or allocation.

## Qualification on 2026-09-09

The release build and 49 single-cell tests in 14 suites passed. All 25 streaming
CLI regression checks passed, including the new excessive-source declaration.
The full Norman product run published and reconstructed successfully, refused
overwrite, and matched all 3,577,278 aggregate values, all 111,445 cell QC records,
feature axes and condition memberships exactly against independent NumPy/h5py
aggregation. Two fresh publishes produced byte-identical reports and receipts;
two reference preparations produced identical arrays, plans and split definitions.

The initial Python comparison was explicitly interrupted after native publish
and verification succeeded: repeated lazy NPZ reads were decompressing arrays
inside its per-cell loop. The checker now loads the reference arrays once; the
complete qualification was rerun successfully. The initial native logs and
interruption history are retained.

Native timing and peak memory are recorded in `performance.json`. The resident
metadata/aggregate/report costs remain significant; these measurements establish
this full 111,445-cell input, not million-cell execution or an end-to-end speed
advantage over scverse. Unseen-perturbation prediction is still unqualified.

## Reproduce

Use the single-cell benchmark Python environment (NumPy and h5py versions are
recorded in `audit.json`) and a fresh release product built from this source.

```sh
curl -fL 'https://zenodo.org/records/13350497/files/NormanWeissman2019_filtered.h5ad?download=1' -o original.h5ad
python Tools/Omics/PerturbationPrediction/Norman/prepare.py --source original.h5ad --out prepared
python Tools/Omics/PerturbationPrediction/Norman/check_native.py --binary /absolute/path/numivivo --source original.h5ad --prepared prepared --out native
python Tools/Omics/PerturbationPrediction/Norman/check_activation.py --prepared prepared --output activation-audit.json
```

Output directories must not exist. The native checker runs the actual product,
reconstructs the bundle, checks overwrite rejection, compares every aggregate
value, all 111,445 cell QC records, all axes and condition membership, and checks
the retained source hash. The independent reference only materializes the small
237-by-33,694 condition matrix, never the cell-by-gene matrix.

The evidence directory retains the reference arrays, mappings, splits, audit,
compressed native report, receipt and validation logs. The 699 MB source is
retrieved using its pinned public URL and checksums rather than committed to Git.
To reconstruct the archived native bundle, restore `report.json.gz` to
`report.json`, use `native-plan.json` as `plan.json`, and download the exact
`original.h5ad`; `receipt.json` belongs to the recorded executable implementation.
Rebuilt executables can publish and verify their own fresh bundles.
