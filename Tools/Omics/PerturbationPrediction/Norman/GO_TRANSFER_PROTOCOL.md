# Gene Ontology transfer: declared development protocol

Declared 2026-09-10 before preparing the complete annotation panel or scoring
this model. Earlier co-response and control-correlation failures inform this
experiment. Norman is a reused development benchmark, not an untouched study.

Use all 105 held-target single-condition folds and all 33,694 original features
from the byte-pinned reference. Exclude every condition involving the held target
from fitting; only controls and other single perturbations are training data.
Do not use paired responses, held responses, post-treatment descriptors or
outer-fold scores to choose annotations or model settings.

## Independent target information

Retrieve GO annotations from MyGene.info's documented gene annotation API, using
the exact original Ensembl feature ID for each uniquely source-matched target.
Retain raw requests/responses, metadata build version, retrieval time, hashes,
Entrez/Ensembl identities, symbol changes, GO namespace, evidence, qualifiers and
PubMed references. Require human taxon 9606 and an exact returned Ensembl match.
Missing, ambiguous or unmatched identities remain unsupported; no alias guessing.
The three source-symbol gaps from earlier experiments remain explicit.

Use distinct directly annotated BP/MF/CC GO identifiers. Exclude annotations with
a NOT qualifier and the three ontology root terms. Retain all evidence classes;
do not infer ancestry, causality, tissue specificity or annotation dates that the
response does not provide. These current annotations can reflect knowledge
published after Norman. Their use does not establish temporal independence or
absence of indirect study overlap. Report PubMed overlap with the Norman paper
when its identity can be verified, without tuning the term set after scoring.

The full target kernel is Jaccard set similarity, intersection over union.
No RNA values enter it. Record coverage, annotation counts and every target pair.
This is a transparent annotation-kernel experiment, not an execution of GEARS.
The GEARS source motivated GO information, but its Dataverse file could not be
retrieved and is not silently substituted or claimed as the input.

## Models and inner selection

Retain no-change, mean of all other singles, and mean of descriptor-supported
training singles. For supported targets, fit kernel ridge with an unpenalized
intercept to log1p-CPM response vectors. Evaluate fixed lambda=1 and a nested
variant selecting lambda from [0.01, 0.1, 1, 10, 100]. Select by the mean
per-training-target, all-gene leave-one-target-out RMSE after clipping predicted
log-expression to zero. Resolve exact ties toward the larger lambda. The outer
held target cannot enter these inner validation outcomes or fits.

Retain a matched shuffled nested kernel model: rotate the ordered supported
training responses by one position relative to GO descriptors. Its inner
selection is repeated on that shuffled training relation. This is a fixed
negative control, not a permutation significance test. Store all candidate inner
losses, chosen lambdas, query weights and training target order.

Clip negative predicted log-expression at zero, retain clipping counts and
implied CPM sums, and do not reclose CPM. Unsupported descriptor predictions are
explicitly absent/NaN; generic no-change and all-single means still score all
105 targets. Freeze predictions and receipts before the scoring command.

## Verification and evaluation

Compare the kernel/intercept solver against an independent centered-kernel
scikit-learn solve. Check analytical inner leave-one-out predictions against
explicit refits. Changing excluded outcomes must leave training inputs and
predictions unchanged. Repeating preparations/predictions/scoring must reproduce
their scientific arrays and receipts exactly; source retrieval timestamps are
archived once and are not represented as deterministic remote state.

Use the original all-gene and training-selected top-1,000 panels, RMSE, MAE,
correlation and sign metrics. The primary comparison is mean per-target all-gene
RMSE against the generic mean, the matched supported-training mean and the
shuffled model, on identical supported outer targets. Report every target and
panel, worse-than-no-change cases and coverage. Improvement must hold against
all three comparators before claiming useful GO-specific transfer on this
development set. No metric threshold or method changes after scoring.

Dense storage is condition/target by gene and target by target; no dense
cell-by-gene matrix. This experiment does not by itself qualify native prediction,
single-cell distributions, unseen tissue/donor transfer, uncertainty calibration,
genetic interactions or Bayesian/mechanistic coupling.
