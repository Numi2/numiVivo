# Focused target and context transfer — 4 October 2026

The experimental v0.4 baseline remains frozen at NumiVivo `3acba615` and NumiLab `cc83364f`. This increment tests one representation replacement and improves the existing shared specimen workspace. Biological promotion remains false.

## Registered experiment and scope

Five nested development folds were executed at fixed native architecture: three globally held-target folds over 46 targets / 81 source–context–target groups, and two held-context folds (32 shared targets each) for stimulated versus unstimulated Jurkat. Complete targets or contexts, source rows, guides and their eight bags stay together. Source-row disjointness and fold-specific preprocessing were checked. The three studies remain human cellular technical development evidence. These are not 81 unseen genes, independent animals, mouse tissue transfer or cross-study biological validation. The original capture6 validation and exposed GSE168620 reserved result were not reused as independent evaluation.

The ESM2 arm uses frozen `facebook/esm2_t33_650M_UR50D`, revision `08e4846e537177426273712802403f7ba8261b6c`, float32 MPS inference, residue mean pooling and length-weighted chunks above 1022 residues. PCA16 is fitted only on each inner/outer training-target set, padded onto the existing descriptor interface. Sequence knowledge about a held target is intentionally available; held response data never fits that representation. No perturbation-response independence claim is made about pretraining. Protein source records, sequence hashes, model-file hashes and the embedding-owner snapshot are retained. [Official model](https://huggingface.co/facebook/esm2_t33_650M_UR50D) · [official ESM repository](https://github.com/facebookresearch/esm).

Controls are the prior composition/GO descriptors and target-ID model, plus matched-input ridge regression. An unknown target has its ID contribution masked, so the held-target ID arm is explicitly a context-only/missing-target control. It is not a transferable biological target representation. Species, CRISPRi/KO and supported stimulated/unstimulated metadata remain explicit. No new architecture or optimizer sweep was run. Native checkpoints 240/1440/5760 and ridge penalties 0.1/1/10 were registered before fitting. Each selector used inner folds, then refit on outer training. The context selector has only three shared K562-high-MOI targets and a mixed study/modality shift; it cannot select a generally validated model.

Common feature identities come from source metadata; variable features come only from training-context controls. Response scaling uses training bags only. Inner and outer-refit feature selection/PCA/scaling are separate. Wrong-target tests exchange descriptors/IDs within source/context while preserving the receiving control profile and context. Original whole-gene response RMSE, direction deadband 0.1 and selection_score definitions remain unchanged. Panels vary by fold, so compare representations within each fold and do not merge the two tasks.

Before any response fit, the protein audit found that the old DDIT3 accession P0DPQ6 was a 34-residue upstream ORF. The new source-bound manifest uses the main DDIT3 protein P35638 (169 residues) for both coarse and ESM2 arms. `protein-mapping-amendment.json` retains the correction; old v0.4 artifacts are unchanged. [DDIT3 reference](https://www.ncbi.nlm.nih.gov/gene/1649).

## Response error by fold

Lower is better. These are equal-group RMSEs, not uncertainty intervals.

| Fold | No change | Coarse native | ESM2 native | ID native | Coarse ridge | ESM2 ridge | ID ridge |
|---|---:|---:|---:|---:|---:|---:|---:|
| held-target-0 | 0.32658 | 0.32221 | 0.37677 | 0.32509 | 0.32349 | 0.34774 | 0.31965 |
| held-target-1 | 0.33747 | 0.34586 | 0.36627 | 0.33889 | 0.33721 | 0.36005 | 0.33265 |
| held-target-2 | 0.37121 | 0.38042 | 0.40214 | 0.37442 | 0.36692 | 0.38484 | 0.36424 |
| held-context-stimulated | 0.43039 | 0.45340 | 0.47393 | 0.49456 | 0.43834 | 0.45957 | 0.43814 |
| held-context-unstimulated | 0.38441 | 0.44897 | 0.48454 | 0.47985 | 0.40359 | 0.43477 | 0.40222 |

## What transferred—and what did not

ESM2 was worse than the coarse descriptor under both native and ridge models in all five folds. It is not adopted. The held-target native coarse model is better than wrong-target inputs on average, but still worse than no change overall. Correct target information has some numerical effect without earning useful predictive error. The masked-ID ridge control has the lowest held-target aggregate error, and its predictions are invariant to exchanging unknown target IDs. Its small gain therefore cannot be credited to unseen-target biology. Direction accuracy is also appreciable for that context-only control, showing why direction alone is insufficient.

| Task / model | Response RMSE | Wrong-target RMSE | Direction accuracy |
|---|---:|---:|---:|
| unseen-target-familiar-context / coarse native | 0.34992 | 0.36059 | 0.5655 |
| unseen-target-familiar-context / coarse ridge | 0.34291 | 0.34201 | 0.5969 |
| unseen-target-familiar-context / ESM2 native | 0.38217 | 0.38452 | 0.5540 |
| unseen-target-familiar-context / ESM2 ridge | 0.36452 | 0.36443 | 0.5692 |
| unseen-target-familiar-context / target-ID native | 0.34657 | 0.34657 | 0.5752 |
| unseen-target-familiar-context / target-ID ridge | 0.33924 | 0.33924 | 0.6157 |
| known-target-unseen-context / coarse native | 0.45118 | 0.46409 | 0.5257 |
| known-target-unseen-context / coarse ridge | 0.42096 | 0.42317 | 0.5117 |
| known-target-unseen-context / ESM2 native | 0.47924 | 0.49675 | 0.5276 |
| known-target-unseen-context / ESM2 ridge | 0.44717 | 0.44680 | 0.5061 |
| known-target-unseen-context / target-ID native | 0.48721 | 0.51091 | 0.5201 |
| known-target-unseen-context / target-ID ridge | 0.42018 | 0.42128 | 0.5063 |

The preceding task summaries weight source–context–target groups, not biological units. No change is 0.34550 in held-target folds and 0.40740 in held-context folds. Every learned configuration is worse than no change in both context tests. ESM2 ridge wrong-target error is essentially unchanged, or slightly better: these results do not establish target-specific transferable prediction. Most selectors chose 240; the held-target-2 ID model chose 1440. Later checkpoints were evaluated but did not win the registered inner criterion. That observation is retained rather than prompting another sweep.

## Intervention selection

`selection-results.csv` retains every fold × representation × estimator × eligible context, with selected intervention, random expected utility, no-intervention gain, actual fitted training-mean ranking, regret and best observed eligible candidate. Candidate sets include no intervention. Counts and eligibility are in each result JSON. No candidate is labeled a reliable winner. The exact objective is reduction of mean log1p(CPM) of FOS/JUN/JUNB/EGR1, with no spatial-preservation claim. It was registered before training; it is not a new independent objective-validation cohort.

The ESM2 native model has mixed held-target selection outcomes, including negative random/no-change gains in stimulated fold0 and unstimulated fold1. In context transfer it selects DHODH for stimulated Jurkat (observed utility −0.4761, regret 0.8477) and TUBB for unstimulated Jurkat (−0.5158, regret 0.7835). Both are worse than random and no intervention, and tie their target-matched training-mean choice. These failures reject useful context-sensitive intervention selection here. Isolated successful choices remain in the report; none overturns the failed qualification.

The decisive next evidence need is a matched-context perturbation study with repeatable control/target responses and resolved biological units. The present experiment separates representation failure from basic optimizer failure: replacing protein features did not help even the linear baseline, and preserved target inputs did not rescue the context task. It does not establish whether the limiting cause is causal target biology, assay/context shift or measurement variability. No new independent outcomes were opened and no universal winner was selected after examining these outer results. Each fold model/selector is frozen in its prediction seal; the registered spatial default stays unchanged.

## Software qualification and reproducibility

All 15 native outer predictions replay with byte-identical named tensor contents. The initial whole-file hash check failed because safetensors key/payload order differs across serialization; its log is retained. Dtypes, shapes and bytes match exactly, with unchanged sealed file hashes. This is numerical qualification only. Thirteen shared-workspace tests and six live owner-evidence checks passed. The old sealed comparison remains unrevealed and byte-identical.

The updated workspace uses one large specimen by default, optional synchronized layers, visual editable cards, local support, explicit unsupported-intent correction, gene-level model-difference attribution and direct failure inspection. Corrected v0.4 and v0.3 are still explicitly selectable. Browser pointer/keyboard and actual browser-driven candidate drag are engineering checks, not an uncoached researcher study. Researcher usability remains pending by user instruction. Physical touch and Safari remain unqualified; native computer-use service startup failed. Touch handler unit tests do not close that gate. No new spatial-preservation claim is made.

Run reproduction on Apple silicon from the retained native runtime and campaign plans. `replay_transfer.py` verifies the original seals and named-tensor replay without editing paths in sealed plans. `transfer_experiment.py register/prepare/run` and `embed_transfer_targets.py --help` document preparation/training; the campaign registration records exact inputs, source hashes, preprocessing, budget and checkpoints. To reconstruct from original h5ad files, use the source identities/hashes and frozen v0.4 cohort release; do not treat relocation as permission to modify scientific content. Frozen ESM2 weights can be obtained from the pinned official revision; the release retains embeddings and hashes rather than repackaging the 2.6GB public checkpoint.

The ESM2 training track and spatial investigation track remain separate. A cellular technical-transfer experiment does not register a new spatial predictor.
