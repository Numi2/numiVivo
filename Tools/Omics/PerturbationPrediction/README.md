# Donor-held-out perturbation response baselines

The independent [Adamson UPR study](Adamson/GEO_RESTORATION.md) now has original
GEO identities restored for all 65,337 cells, correcting 3,634 assignments affected
by its deposited barcode join. Corrected native counts match the independent
reference, and fixed GO descriptors are captured. Its control identities and
prediction scores remain unresolved. The [AlphaGenome Atlas note](ALPHAGENOME_ATLAS.md) records
the later regulatory-variant integration opportunity and its current boundaries.

This is a deterministic external benchmark for the [native prediction
implementation](NATIVE.md). It uses full, source-verified pseudobulk counts from Kang B cells
(2,651 cells, 15,706 genes, eight paired donors) and Hagai fibroblasts (13,863 cells,
22,048 genes, three paired donors). It does not train across species. Every donor
is held out once. Its control profile is supplied, and its treated profile is
used only for evaluation.

The [predeclared protocol](PROTOCOL.md) fixes normalization, training-only context
features and response panels, and four baseline models. All source genes are
retained in the normalization denominator and all genes receive predictions.
No-change, shared mean response, shared median response and control-context ridge
are compared. Ridge is deliberately fixed at alpha=1, with no test-donor tuning.
The treatment is known from other donors; these are not unseen-perturbation tests.

## Results

Mean per-donor response RMSE in natural-log(1+CPM) units; lower is better:

| Study/panel | No change | Mean response | Median response | Context ridge |
|---|---:|---:|---:|---:|
| Kang, all genes | 1.158536 | 1.100834 | 1.101989 | 1.080626 |
| Kang, training top 200 response genes | 2.846039 | 1.364420 | 1.468104 | 1.292907 |
| Hagai, all genes | 0.572022 | 0.247536 | 0.247536 | 0.239330 |
| Hagai, training top 200 response genes | 3.281963 | 0.812301 | 0.812301 | 0.773312 |

The top 200 are selected independently in each training fold and can differ
between folds. Full-gene metrics prevent presenting that selected panel as a
whole-transcriptome result. Ridge improves average RMSE but is not uniformly
better: for mouse1's top-response panel it worsens mean-response RMSE from
1.084105 to 1.132784. On the pre-existing Hagai marker panel, average RMSE worsens
from 0.617506 to 0.625617. With only two training mice per fold, the mean and median
baselines are identical. No fitting warnings occurred.

Saved metrics also report MAE, Pearson correlation when defined, sign agreement
and explained response sum of squares relative to zero change. Expression-space
metrics are separate because baseline gene expression can dominate correlation.
Predictions are clipped to nonnegative logCPM; both the original estimated change
and the applied change are retained. These are compositional pseudobulk expression
predictions, not absolute molecule counts or individual-cell distributions.
The point estimates are not reclosed to a CPM composition: context ridge implies
CPM sums of 1.022-1.096 million in Kang and 0.938-1.541 million in Hagai. The added
composition diagnostic records this deviation; no post-hoc normalization changes
the predeclared predictions or scores. These outputs must not be rounded or
relabelled as calibrated count libraries.

## Reproduction and evidence

```sh
python verify_inputs.py --input /prepared/kang --out kang-input-check.json
python run.py --input /prepared/kang --control ctrl --treated stim --study Kang-B-cells --out /new/kang
python run.py --input /prepared/hagai --control unstimulated --treated LPS6 --study Hagai-fibroblasts --out /new/hagai
# Repeat into new directories, then compare arrays and metadata exactly:
python check_replay.py /new/kang /new/kang-repeat --out kang-replay.json
```

The input directories contain `counts.tsv`, `samples.tsv` and `input.json` from
the existing independent Bioconductor comparison preparation. Preserved count
TSVs are gzip-compressed in this evidence tree; decompress before running.
`verify_inputs.py` needs the original H5AD and its mapping at the locations recorded
in `input.json`, or those locations updated explicitly after restoring the same
hash-verified source. It aggregates sparse cells to a small donor-by-gene matrix
and verifies all raw counts, axes, library totals, donors and conditions.

Both studies passed raw H5AD reaggregation and exact repeated-run comparison.
Model construction accepts only training matrices. Changing the sealed donor's
treated counts leaves fitted arrays and predictions unchanged in all 11 folds.
Frozen model reload reproduces predictions. scikit-learn dual coefficients and
predictions agree with the independent NumPy solve at 1e-10 tolerances. NPZ arrays
use Unicode IDs and load with `allow_pickle=False`.

Evidence includes all input count tables and metadata, per-fold model parameters,
feature panels, donor membership, per-gene predictions and outcomes, full metrics,
package versions, source hashes, source reaggregation and replay checks. The
original single-cell source identity is inherited from the existing experimental
benchmark artifacts. The original external benchmark commit introduced no
new native executable; the native owner is qualified separately.

## Native integration boundary

The native contract now binds a training cohort and perturbation identity,
complete feature universe/namespace, assay units, training-only transforms, model
parameters and allowed query inputs. Treated outcomes remain outside prediction.
It reuses the Omics source/aggregation/provenance owners and preserves all four
baselines. See the native qualification for the exact supported scope.

The existing `Calibration/VivoPosteriorModel.swift` and
`VivoTemperedPosteriorSampler.swift` already provide likelihood fingerprints,
bounded priors, explicit evaluation budgets, checkpoint identity and batched
likelihood evaluation. An omics adapter can supply an authoritative likelihood;
it must preserve failed/nonconverged populations and independently validate
predictive intervals. The existing negative-binomial owner supplies count-model
arithmetic, but neither this ridge benchmark nor the current bounded-parameter
sampler is automatically a hierarchical whole-transcriptome predictor. Learned
omics-to-mechanistic parameter mappings and physical units need separate models
and calibration. No Bayesian or mechanistic prediction is claimed here.

Bayesian/mechanistic prediction, unseen perturbation identities, unseen cell/tissue
contexts, heterogeneous single-cell responses, uncertainty calibration and larger
independent donor cohorts remain open.

## Rationale and primary sources

Strong shared-response and linear baselines are necessary because expression
correlation can hide weak perturbation-specific prediction; see the primary
[linear-baseline comparison](https://www.nature.com/articles/s41592-025-02772-6)
and [Systema evaluation](https://www.nature.com/articles/s41587-025-02777-8).
[scGen](https://www.nature.com/articles/s41592-019-0494-8) demonstrates a different
latent-space perturbation-prediction approach; it is not run or outperformed by
this benchmark. The paired donor studies are the same source-qualified Kang and
Hagai inputs used in [the Bioconductor comparison](../Bioconductor/README.md).
