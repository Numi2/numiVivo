# Adamson provisional role sensitivity

**Status: completed as a provisional sensitivity experiment; not an
authoritative biological qualification.** The primary paper and Table S1
identify the UPR control set, but the restored GEO labels `63(mod)_pBA580` and
`Gal4-4(mod)_pBA582` are not mapped to `NegCtrl-2` or `NegCtrl-3`. This run
therefore pools those two deposited labels under an explicit provisional
control condition and keeps the role gate closed.

## What ran

The run used the restored selected Adamson cohort: 50,440 cells, 94 deposited
guide groups, 32,738 genes and 82 candidate target prefixes with positive
selected counts. All deposited guide rows present for a target prefix were
aggregated together before fitting. Eighty prefixes had captured GO terms;
`IER3IP1` and `TIMM23` were retained as no-data descriptor cases and were not
included in the primary score.

Each of the 82 targets was held out in turn. The native target-kernel commands
fit on the provisional control plus the other target aggregates, verified the
model, predicted the held descriptor and verified the prediction. Every
prediction was frozen with its query and receipt hashes before the separate
score pass read the held aggregate. The algorithm was unchanged: log1p-CPM,
direct-term Jaccard similarity, ridge lambda 1 with an unpenalized intercept,
zero clipping and no CPM reclosure.

## Conditional result

The primary endpoint is mean per-target all-gene log-response RMSE over the 80
supported held targets:

| Method | Mean RMSE | Held targets |
| --- | ---: | ---: |
| GO-term kernel | **0.1157221444** | 80 |
| All-training-single mean | 0.1221013658 | 80 |
| Supported-training-single mean | 0.1221668265 | 80 |
| Matched supported-response shuffle | 0.1229151825 | 80 |
| No change | 0.1300771830 | 80 |

Under this declared pooled-control sensitivity, the GO-term kernel is lower
than all three predeclared learned baselines in mean RMSE. It beats the
all-training-single mean in 64/80 targets, the supported-single mean in 64/80,
the matched shuffle in 63/80 and no change in 50/80. These are conditional
expression-response scores for one pooled K562 source; they do not identify a
control construct or establish transfer to another donor, tissue, disease,
clinical endpoint or phenotype.

## Evidence and limits

The immutable remote evidence root is
`/Users/n/numivivo-adamson-provisional-role-sensitivity-20260914-retry` on the
physical Mac mini. Its compact receipt and score identities are:

| Artifact | SHA-256 |
| --- | --- |
| Native selected-cohort report | `fd23394c260a1be7daecfc692a444cdef6b858b8bc7f55b89192bf559fae9a78` |
| Selected-cohort receipt | `e8a8ebf59cafc09f2313059825eec3c71f6e8c33e6945cc6f1670eff3761a126` |
| GO coverage | `1dc7b844c101509be36905fd8aeae4478fae401482ef92eea6a529efea51e13e` |
| GO annotations (`.json.gz`) | `f867977cb4fa170df804ff8251421f2be3709df4625df84a8ee8c47dbd09726f` |
| Native executable | `93986de58266a8cbc6f43f51537907ad1e39ba8a09a0bdbbfdfcdf95f9b72469` |
| Provisional training JSON | `9ef74639c36ba2ed7bd8e2c881944d3bb933b41df1b160ac08351cef64c57a17` |
| Driver | `6672818d2e742b0801d9ba83b7aa8298d5fd13927f8b3d14f44c18fa8e4a3a5c` |
| Provisional receipt | `7752f65e8182d0ff7673c781a7131ce4e3a40a8c06967ce6f3acb5611f674b0a` |
| Score summary | `31253574aee5c9575a010ec5f31b3fc6c7ce9fdb167facdf143eeb214da6b184` |

All 328 native commands returned zero. The 82 retained compressed reports,
queries and prediction receipts match the frozen hashes, and no native staging
directory remains. The first full run is retained separately as a diagnostic;
it had a reporting-only mismatch between an uncompressed report hash and its
retained gzip file, then was rerun with the corrected retained-artifact hash.

The result does not change `controlsVerified`, `experimentalTargetsVerified`,
or the frozen Adamson qualification status. To promote this from sensitivity
evidence, obtain an authoritative deposited-label-to-sequence mapping and
reconcile the 94 selected groups with the paper's 93-guide summary, then rerun
the frozen protocol without using these scores for tuning.

## Reproduce

The standard-library driver is
[`run_provisional_role_sensitivity.py`](../Tools/Omics/PerturbationPrediction/Adamson/run_provisional_role_sensitivity.py).
On the Mac mini, with the retained selected cohort and descriptor files:

```sh
python3 run_provisional_role_sensitivity.py \
  --binary /Users/n/numivivo-adamson-target-kernel-20260910/numivivo \
  --cohort-report /Users/n/numivivo-adamson-target-kernel-20260910/cohort/native/report.json \
  --cohort-receipt /Users/n/numivivo-adamson-target-kernel-20260910/cohort/native/receipt.json \
  --coverage /tmp/adamson-coverage.json \
  --annotations /tmp/adamson-annotations.json.gz \
  --out /Users/n/numivivo-adamson-provisional-role-sensitivity-<new-run>
```

The output is explicitly labelled provisional even when every numerical and
reconstruction check passes. It establishes a conditional RNA-response
estimate from available counts and annotations; it does not establish a
general biological-outcome predictor.
