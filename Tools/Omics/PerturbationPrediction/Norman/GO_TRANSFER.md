# GO-informed unseen-target transfer: modest fixed-model gain

A fixed Gene Ontology kernel improves the average held-target response error over
both mean-response baselines and its matched shuffled control on the full Norman
development benchmark. The nested model fails the primary all-gene comparison.
Both models, their choices and every failed target remain recorded. This is a
candidate for a native implementation and independent-study validation, not a
qualified general unseen-target predictor.

The [primary protocol](GO_TRANSFER_PROTOCOL.md) was frozen before preparing the
complete annotation panel or scoring this model. The matched fixed-lambda shuffle
was a separately declared [post-result audit](GO_FIXED_CONTROL_PROTOCOL.md), added
because the primary shuffled model used nested selection. It was not part of the
original predeclared comparison. Neither audit changed any prediction parameter.

## Source identity and coverage

The byte-pinned Norman condition reference remains
`ec22f61196cf51f7c1a62c681e4076bd9ea269429f1423dfd65e21c11a7ce8bf`.
All 105 single-target folds and all 33,694 source features remain. Each fold removes
every single or paired condition involving its target; fitting uses only controls
and other singles. The fitter never receives the held response. The orchestrator
reads the reference to select these arrays; this is function-level separation,
not an operating-system sandbox.

The [GEARS paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11180609/) motivates
using GO information for perturbation identity. Its pinned source revision
`f374e43e197b295016d80395d7a54ddb81cc6769` points to Dataverse file 6153417, but
that download returned HTTP 403. The failed request is retained; neither GEARS
weights nor its dataset were used.

Instead, the [MyGene.info annotation API](https://docs.mygene.info/en/latest/doc/annotation_service.html)
returned annotations for exact original Ensembl gene IDs. Both metadata snapshots
report build `20260906`; the Entrez source reports `20260904`, and UniProt reports
`20260610`. The complete raw responses, URLs, retrieval times, metadata, hashes,
GO records and original/current gene identities are retained. Human taxon 9606
and the exact returned Ensembl ID are checked before use.

There are 101 supported targets and 1,259 distinct GO terms. C19orf26, C3orf72
and KIAA1804 retain their unresolved source-symbol status. IER5L maps correctly
but has no usable GO annotations. No aliases are guessed and no expression
threshold removes weakly expressed targets. Generic baselines still cover 105.

The kernel is Jaccard similarity of distinct directly annotated BP/MF/CC terms.
NOT-qualified annotations and the three ontology roots are excluded; ancestry is
not imputed. All evidence classes remain, with qualifiers and citations archived.
Independent binary-matrix intersections reproduce every set-based kernel value
exactly. The minimum eigenvalue is 0.293838, so this panel is positive definite.

None of the 4,034 retained annotation records directly cites the
[Norman paper, PMID 31395745](https://pubmed.ncbi.nlm.nih.gov/31395745/).
That check does not establish temporal independence or exclude indirect study
overlap. These are current annotations, reused development data and pooled K562
CRISPRa means, without independent biological replication.

## Model and verification

Kernel ridge fits all-gene log1p-CPM response vectors with an unpenalized intercept.
The fixed model uses lambda=1. The nested model chooses among 0.01, 0.1, 1, 10,
and 100 by all-gene leave-one-training-target-out RMSE, with larger-lambda tie
resolution. Clipping is applied before inner scoring. No outer outcome enters
selection. Negative predicted log-expression is clipped to zero, with no CPM
reclosure; clipping counts and implied CPM sums are retained.

For K and lambda, let C=(K+lambda I)^-1, u=C1, and
A=C-uuᵀ/(1ᵀu). Query weights are kᵀA+uᵀ/(1ᵀu), which sum to one. Inner
leave-one-out predictions use I-A/diag(A), with row-wise diagonal division.
An independent centered-kernel scikit-learn solver checks final predictions.

Two preparations and two complete prediction runs are byte-identical, including
receipts. Independent prediction error is at most 2.00e-15. There are 45 synthetic
explicit leave-one-out refits (maximum error 2.14e-13) and 30 real-data explicit
refits (2.14e-15). Every fold's selected input arrays are invariant to changing
all excluded outcomes. Full mutated refits for AHR, ZNF318, IER5L and KIAA1804
reproduce predictions, selected panels, weights and inner choices exactly; the
last two exercise distinct unsupported cases.

The matched fixed shuffle reuses the exact frozen fixed-model query weights,
rotates only supported training responses, and agrees with an independent
lambda=1 refit within 1.45e-15. Both its complete runs and both scoring runs are
byte-identical. All 105 target identities and all unsupported NaN rows are checked.

## Held-target results

Mean per-target response RMSE on the same 101 supported targets:

| Method | All genes | Training top 1,000 | Worse than no-change, all genes |
|---|---:|---:|---:|
| No change | 0.133391 | 0.122163 | 0 |
| Mean of all other singles | 0.127227 | 0.105720 | 23 |
| Mean of supported training singles | 0.127239 | 0.105737 | 22 |
| GO kernel, fixed lambda=1 | 0.126490 | 0.102073 | 29 |
| GO kernel, nested lambda | 0.127325 | 0.103728 | 29 |
| Shuffled GO, nested lambda | 0.127010 | 0.105097 | 22 |
| Shuffled GO, fixed lambda=1, post-result control | 0.127379 | 0.104376 | 26 |

Fixed GO improves average all-gene RMSE by 0.000737, approximately 0.58%, versus
the all-single mean. It improves 59/101 targets on all genes and 65/101 on the
training top-1,000 panel. Against its matched fixed shuffle it improves 57/101
and 63/101, respectively. The gain is small and nonuniform: 42 targets are worse
than the mean baseline, and 29 are worse than predicting no change. These are
engineering comparisons, not confidence intervals or biological replication.

Nested selection chooses lambda=1 in 90 folds and 10 in 11; its shuffled control
chooses 10 in 97 folds and 100 in four. The nested model loses to both mean
baselines and its matched shuffle on all genes, despite better mean top-panel
error. The primary failure is retained, not replaced with the secondary panel.
A training-only selection rule therefore did not establish better outer-target
transfer in this experiment.

Earlier co-response and control-only experiments used 102 and 97 supported
cohorts. Their published means cannot be directly compared with these 101-target
means as if coverage were identical. The present comparisons use matched targets.
Full per-target RMSE, MAE, response correlations and sign metrics are retained.

The repeated primary prediction run took 23.49 seconds and peak RSS
1,154,908,160 bytes on the local Apple M4. This includes all 105 folds, inner
selection, independent solvers, isolation checks and compression. The data are
condition-by-gene, not dense cell-by-gene. These are observed reference-run costs,
not native speedup, Metal or million-cell qualification.

## Reproduce and remaining work

`go_transfer.py` exposes fetch, prepare, predict and score commands. Use the
archived annotation capture to reproduce the exact prepared kernel; a new live
fetch can reflect a different service build. Paths are explicit through `--help`.
Run preparation and prediction twice, then `check_go_transfer.py`, and invoke
scoring only after prediction receipts exist. `go_fixed_control.py` separately
reproduces the post-result control using the original primary artifacts.

The archive in `evidence/2026-09-10-go-transfer` retains raw annotation sources,
prepared kernels, both predictions, per-target scores, protocols, declarations,
checks and logs. Large prediction archives are split into exact byte chunks to
keep individual Git objects bounded; the chunk map records whole-file and chunk
hashes. Reassembly is byte-for-byte, not a new serialization of the predictions.

The archive verifies 167 logical files in 141 unique compressed objects,
277,465,048 logical bytes and 135,843,124 stored object bytes, plus 12 readable
summaries. Both large prediction files reconstruct to their original hashes.
From this directory, run:

```sh
python archive_go_transfer.py --verify --out evidence/2026-09-10-go-transfer
```

For restoration, decompress each `manifest.json` file entry's `object` to its
logical `path`, then concatenate each `largeArtifacts` entry's ordered `parts`
to its `path`. Verify the recorded byte count and SHA-256 before loading it.
The original qualified Norman reference stays external with its exact hash.

The useful next implementation candidate is the explicit fixed kernel, while the
nested failure rules out promoting this selection rule. Native reusable model
bundles, independent untouched-study validation, unknown biological contexts,
calibrated uncertainty, single-cell response distributions and Bayesian/mechanistic
coupling remain open. No Swift behavior or native product qualification changed
in this reference experiment.
