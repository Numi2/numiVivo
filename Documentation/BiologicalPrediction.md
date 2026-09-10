# Can NumiVivo predict biological outcomes?

**Yes, within specific gene-expression experiments. Reliable prediction of
general biological outcomes is not established.** The strongest evidence is
prediction of a known perturbation in an unseen donor from that donor's observed
control profile. Native models also predict combinations of previously observed
targets, and a GO-based prototype predicts some unseen targets with a small,
nonuniform advantage over simple baselines.

This assessment reviews the expression-prediction evidence and subsequent
full-cohort scoring, clustering and decoder-calibration results through
`8f4a6f9d`, with an experimental-role audit on 2026-09-10.
It distinguishes measured held-out expression outcomes from numerical
reconstruction, integration diagnostics and conditional molecular simulations.
Each linked experiment retains its actual source, executable and platform
identities; this review does not requalify historical receipts under a new build.

## What information is sufficient for the implemented predictors?

| Question | Information available before prediction | Output and present evidence |
| --- | --- | --- |
| How will a new donor respond to a known treatment? | Paired control/treated count profiles from training donors; the new donor's control counts; matching treatment, organism, population and feature identities. | Native no-change, mean, median and context-ridge estimates of population-average treated log1p(CPM). Complete Kang, Hagai and HIRISA held-out donor experiments show conditional predictive signal, with failures retained. |
| What happens when two observed targets are perturbed together? | Control and single-target counts for both constituents in the same experimental context; the pair's target identities. | Six native composition baselines evaluated on all 131 held-out Norman pairs. This tests expression composition; it does not identify genetic interactions. |
| What happens when an unobserved target is perturbed? | Control and other single-target counts in the same context; independently supplied target identities and GO terms for training and query targets. | Native target-kernel prototype, evaluated on 105 Norman folds with 101 supported descriptors. Its small average gain does not establish reliable unseen-target transfer. |
| What will happen in a new tissue, species, disease state or patient? | Would require a validated transfer model and measured outcomes in that destination context. | No qualifying result in this evidence set. Existing point estimates cannot be promoted to those outcomes. |
| Can a DNA variant predict a cellular or tissue phenotype? | Would require assembly/allele-resolved regulatory evidence and validated connections through RNA, proteins, mechanisms and phenotype. | [AlphaGenome integration work](AlphaGenomeAtlas.md) provides a direction and separate evidence interface; a validated end-to-end phenotype predictor is not established. |

The expression outputs are compositional pseudobulk point estimates, not absolute
molecule counts or individual-cell response distributions. Negative predicted
log-expression is clipped to zero; implied CPM totals are retained rather than
silently renormalized. These models supply neither calibrated predictive intervals
nor an experimentally validated RNA-to-phenotype mapping.

## What the measured outcomes show

### Known treatment, held-out donor

The [HIRISA experiment](../Tools/Omics/Benchmarks/HIRISA/PREDICTION_RESULTS.md)
uses all 79 frozen folds across sixteen population/treatment contrasts, with
18,082 genes. Each fit has four training donor pairs, except Bcell IFNg, which
has three. The held-out donor's treated counts enter scoring only; training
selection, scaling and alpha-one ridge fitting use training data. Predictions
were frozen before scoring. Native reconstruction and independent NumPy checks
pass for every fold.

Mean, median and ridge each improve full-gene contrast RMSE over no-change in
**14 of 16 contrasts**. Ridge improves on the training-mean response in only
**4 of 16**. Monocyte IFN-L1 and NK IFNg are worse than no-change under every
learned baseline. Thus, much of the measured predictability comes from a shared
treatment response; the evidence does not justify presenting donor-context ridge
as consistently superior or selecting a production winner after seeing results.

Examples below are equally weighted held-out donor means, in natural-log(1+CPM)
units; lower response RMSE is better. The linked report retains all sixteen
contrasts and both feature families, including weak and negative results.

| Population / treatment | No change | Mean response | Context ridge |
| --- | ---: | ---: | ---: |
| Bcell / IFNa | 0.274165 | 0.097723 | 0.098362 |
| Monocyte / IFNg | 0.564236 | 0.282629 | 0.261827 |
| Monocyte / IFN-L1 | 0.181475 | 0.188058 | 0.193165 |
| NK / IFNg | 0.085704 | 0.091158 | 0.094143 |

The [Kang and Hagai experiments](../Tools/Omics/PerturbationPrediction/README.md)
add eight and three held-out donor folds respectively. Full-gene mean response
RMSE for no-change / mean / ridge is 1.158536 / 1.100834 / 1.080626 for Kang
B cells and 0.572022 / 0.247536 / 0.239330 for Hagai fibroblasts. These are separate
within-study experiments, not training across species. Weak individual folds
and marker-panel regressions remain reported. Millions of measured cells do
not turn a handful of donors into millions of independent biological replicates.

### Held-out combinations of observed targets

The [Norman composition experiment](../Tools/Omics/PerturbationPrediction/Norman/COMBINATIONS.md)
withholds all 131 paired conditions together and trains on control plus 105
single-target conditions. Across all 33,694 genes, mean constituent response
has mean RMSE **0.162827**, compared with **0.192617** for no-change and
**0.180105** for the mean-single baseline. It is worse than no-change for six
pairs. Full log-additivity is worse than no-change for sixty pairs; its stronger
selected-panel result does not erase that failure.

The [native owner](../Tools/Omics/PerturbationPrediction/Norman/NATIVE_COMPOSITION.md)
reproduces all 786 pair/method vectors. This same-cell-line pooled experiment
does not provide independent biological-replicate uncertainty or qualify
unseen constituent targets, donor transfer or causal interaction effects.

### Unseen targets: implemented, with weak development evidence

The [native GO target kernel](../Tools/Omics/PerturbationPrediction/Norman/NATIVE_TARGET_KERNEL.md)
withholds each of 105 targets, using control and the other singles. Four targets
have unavailable descriptors and retain generic baselines only. On the same
101 supported targets, all-gene mean RMSE is **0.126490**, versus **0.127227**
for the all-training-single mean and **0.127379** for the matched fixed shuffle.
The approximately **0.58%** mean improvement is small: the kernel is worse than
the mean for **42/101** targets and worse than no-change for **29/101**.

These results use repeatedly inspected Norman development data and current GO
annotations without demonstrated temporal independence. The nested
regularization experiment fails the primary all-gene comparison; earlier
co-response and control-descriptor models also fail to beat their simple
baseline. Numerical agreement on 517 native prediction vectors verifies the
implementation, not generalization to a new experiment.

The independent [Adamson preparation](../Tools/Omics/PerturbationPrediction/Adamson/COHORT.md)
has restored deposited cell-to-guide identities and verified its original-author
assignment cohort. The [experimental-role audit](../Tools/Omics/PerturbationPrediction/Adamson/EXPERIMENTAL_ROLES.md)
reconstructs all 50,440 selected cells from original GEO records, but controls and
guide-to-target roles remain unverified. Its 94 selected guide groups also need
reconciliation with the paper's stated 93-guide experiment. No independent
prediction scores are claimed. Primary experimental roles and the complete
roster must be resolved before this fixed algorithm is fitted in Adamson.

## What the million-cell work establishes

Complete HIRISA ingestion, count-based DE, PCA, neighbor graph and native seed-7
integration have published operational and numerical evidence. The
[integration result](../Tools/Omics/Benchmarks/HIRISA/FULL_INTEGRATION_RESULTS.md)
and three full Harmony references pass the frozen coarse response-preservation
margins. The subsequent [measured-program diagnostic](../Tools/Omics/Benchmarks/HIRISA/INTEGRATION_PROGRAMS.md)
finds three control-sensitive preservation failures in native integration and
each of three Harmony runs. Eighteen of 32 comparisons lack sufficient erasure
control sensitivity; some baseline decoders have negative within-library skill.
Every fold was independently checked with per-cell weighted SVD regressions.
A separately declared [within-library fitting follow-up](../Tools/Omics/Benchmarks/HIRISA/PROGRAM_CALIBRATION.md)
then increases sensitivity to 28/32 comparisons and meets the same loss margins
in all four candidates. Direct per-cell SVD checks verify every fit. This shows
that the earlier losses depend on the decoder objective; it is development
evidence after known outcomes, not independent biological replication. Four
controls remain insufficient, and complete preservation remains unqualified. Rare-cell and native multi-seed preservation
remain open. Full-cohort [clustering publication, native replay and independent
partition checks](../Tools/Omics/Benchmarks/HIRISA/FULL_CLUSTERING_RESULTS.md) now
pass for all original cells; graph communities are not validated cell types.

Integration classifiers use full-cohort preprocessing and assess retained
information. They are not prospective prediction of an unseen donor's treated
outcome. Donor-response prediction above uses the separate frozen count-based
pipeline and does not inherit leakage or qualification from those classifiers.
Lower donor-associated variance also does not isolate technical batch removal.

## Next evidence needed

1. Retain the completed full-cohort clustering evidence and the separate original
   mapped-reader baseline under their actual identities. Extend biological
   validation without treating graph partitions as authoritative annotations.
2. Validate the revised within-library diagnostic on independent contexts and
   measured endpoints, retaining four insufficient controls and the original
   decoder/Kang NK-recall failures. Complete rare-cell and neighborhood checks.
   The [compact native program artifact](../Tools/Omics/Benchmarks/HIRISA/NATIVE_PROGRAM_RESULTS.md)
   now passes full-cohort publication, replay and independent comparison.
3. Resolve Adamson controls and its 94-versus-93 guide roster from primary records, then execute the
   frozen independent-study target-prediction protocol with coverage, all
   failures and matched simple/shuffled baselines. Do not tune it on test scores.
4. Execute the [frozen HIRISA preparation-transfer split](../Tools/Omics/Benchmarks/HIRISA/CONTEXT_TRANSFER.md):
   sixty cross-preparation folds and sixty matched within-preparation references.
   Membership and independent sparse counts pass; native aggregation and
   predictions remain pending. All source cells and excluded strata are retained.
5. Establish predictive interval coverage and useful improvements over simple
   baselines on independent biological replicates before promoting a predictor.
   Connecting expression to a measured phenotype requires its own model,
   quantitative units and held-out outcome experiment.

The [single-cell roadmap](SingleCellInteroperability.md) retains the complete
twelve-part objective. GPU acceleration and broader cross-scale biology follow
stable algorithms and biological acceptance; source availability alone cannot
close those gates.

## Use and verify

Use [native donor-response commands](../Tools/Omics/PerturbationPrediction/NATIVE.md),
[composition commands](../Tools/Omics/PerturbationPrediction/Norman/NATIVE_COMPOSITION.md)
or [target-kernel commands](../Tools/Omics/PerturbationPrediction/Norman/NATIVE_TARGET_KERNEL.md)
according to the available inputs in the first table. Keep the complete feature
universe, explicit control/perturbation identities and output provenance.

For this assessment, the HIRISA contrast counts were recomputed from all archived
fold metrics and checked against the stored summaries. Its score JSON SHA256 is
`8d49cd4eb3ce2ced43f006cec931d848903251cfff0d496621004cc3637ce44f`.
All 429 HIRISA prediction archive members and 1,292 logical native target-kernel
files were rechecked for stored and decoded identities during this review.
These are evidence-integrity checks, not newly fitted experiments.

```sh
python Tools/Omics/Benchmarks/HIRISA/verify_archive.py \
  Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-10-prediction
python Tools/Omics/PerturbationPrediction/Norman/archive_native_target_kernel.py \
  --verify --out Tools/Omics/PerturbationPrediction/Norman/evidence/2026-09-10-native-target-kernel
```
