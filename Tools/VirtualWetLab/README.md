# Virtual Wet Lab: measured cell-response assay

This is the NumiVivo owner adapter for the Virtual Wet Lab workspace in NumiLab.
It composes the existing native `singlecell-perturbation-batch` predictor and
its raw-count replay verifier. Python handles experiment authoring and held-out
scoring; it implements no alternative prediction model.

## Prepare and run

Use a native CLI exposing `singlecell-h5ad-pseudobulk`,
`singlecell-perturbation-batch`, and `singlecell-perturbation-batch-verify`.
Set `NUMIVIVO_HDF5_LIBRARY` to the native HDF5 dynamic library. Preparation needs
Python with anndata, numpy, scipy and pandas; the workspace adapter itself uses
only the Python standard library.

```sh
curl -fL https://ndownloader.figshare.com/files/34464122 -o kang_2018.h5ad
export NUMIVIVO_HDF5_LIBRARY=/absolute/path/libhdf5.dylib
python Tools/VirtualWetLab/prepare_kang.py kang_2018.h5ad \
  --binary /absolute/path/numivivo --output /absolute/path/kang
python Tools/VirtualWetLab/wetlab.py catalog /absolute/path/kang/assay.json
numi wet-lab --vivo-root "$PWD" --binary /absolute/path/numivivo \
  --assay /absolute/path/kang/assay.json --workspace /absolute/path/experiments
```

The preparation accepts only the registered public H5AD SHA256
`e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830`.
It retains all 15,706 genes and all 2,651 source-labelled B cells from eight
paired donors, creates an annotation-only measured-count projection, and retains
the original source, row mapping and preparation hash. No genes are selected by
held-out outcomes. The study reports six-hour IFN-beta stimulation; no dose or
time interpolation is available. Primary-library accession mapping is unresolved
and batch is unreported. [Source study](https://doi.org/10.1038/nbt.4042).

For command-line use:

```sh
python Tools/VirtualWetLab/wetlab.py predict /absolute/path/kang/assay.json \
  --binary /absolute/path/numivivo --workspace /absolute/path/experiments \
  --donor patient_101 --hours 6 --intervention IFN-beta
python Tools/VirtualWetLab/wetlab.py inspect /absolute/path/experiments/RUN_ID
python Tools/VirtualWetLab/wetlab.py reveal /absolute/path/experiments/RUN_ID \
  --binary /absolute/path/numivivo
python Tools/VirtualWetLab/wetlab.py verify /absolute/path/experiments/RUN_ID \
  --binary /absolute/path/numivivo
```

Prediction writes registration before fitting, uses seven other donors, and
supplies only the query donor's untreated aggregate to the native predictor.
It seals all native artifacts before scoring. Reveal first replays the native
bundle, then computes all-gene response RMSE/MAE against the held-out treated
aggregate in natural-log(1+CPM). Context ridge must strictly beat no-change and
training mean for the displayed per-donor RMSE verdict. Median is also retained.
There is no statistical-significance or biological-qualification claim.

The native source snapshot contains the later observations; this is computational
holdout and UI blinding, not access-controlled escrow. Every failure retains its
logs. Repeat runs are deterministic replays, not independent replicates. Hash
receipts detect accidental changes and the verifier reconstructs predictions;
they are not third-party signatures. Retain the full experiment directory and
exact CLI, HDF5 library and adapter for replay. JSON downloaded from NumiLab is a
summary, not a substitute for raw source artifacts.

## Acceptance

```sh
python -m unittest discover -s Tools/VirtualWetLab -p 'test_*.py' -v
python Tools/VirtualWetLab/check_real.py --assay /absolute/path/kang/assay.json \
  --binary /absolute/path/numivivo --output /absolute/path/new-acceptance
```

The bounded acceptance freezes all eight donors and a 1e-9 numerical tolerance
before scoring. It compares native pseudobulks against original measured counts,
fits an independent NumPy ridge solve, checks all four model vectors and RMSE,
and exercises native replay, unsupported selections and prediction tampering.
Unit regressions cover duplicate arms, aliased biological units, shared cells,
known-value scoring and source tampering.

[October 3 checks](evidence/2026-10-03/checks.json): all eight donor workflows
passed, with 15,706 genes per donor. Context ridge improved RMSE over training
mean by 0.074–4.065% and beat no change for each donor. These are modest results
on a previously inspected development cohort. They establish neither fresh
biological validation nor phenotype, survival, causal or clinical prediction.
The acceptance reuses the existing predictor; no model parameters were tuned
using these scores. See [BiologicalPrediction.md](../../Documentation/BiologicalPrediction.md)
for the broader suite's failures and evidence limits.

## New assay families

The v1 adapter explicitly admits only `family: cell-response`. Add molecular or
spatial-tissue adapters when their native owners provide supported conditions,
controls, measured endpoint, unit-aware scoring and replay. Do not route them
through this RNA scorer. NumiLab's `docs/VIRTUAL_WET_LAB.md` specifies that boundary.
