# Can NumiVivo predict biological outcomes?

**Available information supports conditional RNA-response estimates with limited,
context-dependent predictive value. Reliable prediction of general biological
outcomes is not established.** The strongest evidence is
prediction of a known perturbation in an unseen donor from that donor's observed
control profile. Native models also predict combinations of previously observed
targets, and a GO-based prototype predicts some unseen targets with a small,
nonuniform advantage over simple baselines.

The practical decision depends on the requested outcome and what was observed
before it. For a known treatment, paired training-donor RNA measurements and
the query donor's untreated RNA can support a conditional expression estimate.
For an unseen target, gene identity and annotations are additional inputs;
the fixed GO model now has a modest separately collected UPR validation result,
but reliable generalization across targets and biological contexts remains unproven. DNA or RNA
information alone does not supply a validated mapping to survival, tissue
function, disease progression or treatment benefit. Those outcomes require
their own measured endpoints and held-out validation.

This assessment reviews the expression-prediction evidence and subsequent
full-cohort scoring, clustering, decoder calibration, preparation transfer and
annotation retention through 2026-09-11, with the retained Adamson experimental-role audit.
It distinguishes measured held-out expression outcomes from numerical
reconstruction, integration diagnostics and conditional molecular simulations.
Each linked experiment retains its actual source, executable and platform
identities; this review does not requalify historical receipts under a new build.

## Current completion and admission gates

This includes the published `71ce4fff` implementation evidence and the subsequent
complete GSE181897 test using its separately identified native executable.
A completed execution check and a passed biological
comparison answer different questions. An unscored cohort is neither a prediction
success nor a failure.

| Requirement | Current evidence | What remains before a stronger claim |
| --- | --- | --- |
| Preserve and analyze real counts | Full Parse ingestion and source replay pass; paired native DE and all six edgeR/limma/DESeq2 comparisons complete. Analysis peaks at 202.3 MiB with unchanged report bytes. | Other pipeline stages retain separate memory bounds; numerical agreement does not establish biological accuracy or false-discovery calibration. |
| Preserve biology through integration | Complete HIRISA execution and coarse response checks pass. Native annotation retention fails 37/146 supported, sensitive comparisons, including 10/29 rare comparisons. | Resolve preservation losses and validate independently; do not promote the protected-stratum candidate, which still fails 37/146. |
| Predict an unseen target | Replogle's fixed GO model passes its aggregate primary comparison in all five technical groups. It loses to the mean in 37/150 target/group folds. | Independent biological contexts, target selection and replication; technical groups do not supply these. |
| Transfer RNA response across studies | HIRISA-trained mean improves GSE181897 RMSE by 5.62%, passing the frozen 5% target across 62 query donors; Kang-trained mean fails. GSE226572 and Kang–HIRISA failures stand. | Establish reproducible utility across contexts and training origins; one origin's pass does not erase another's failure. |
| Quantify predictive uncertainty | Mean-response intervals are implemented and assessed. In GSE181897, HIRISA-trained nominal 95% treated-expression coverage averages only 35.54%; Kang coverage averages 92.39% with much wider intervals. | Independent calibration and useful width; both missing features and donor-level undercoverage remain explicit. |
| Predict tissue, disease or treatment outcomes | No validated RNA/variant-to-endpoint chain is established by these experiments. | Explicit measured endpoints, models linking the quantities and held-out outcome validation. |

The GSE181897 input contract is now resolved and its frozen comparison is complete.
The remaining prepared cohorts and the completed external test are:

| Prepared cohort | Ready | Required before fitting or scoring |
| --- | --- | --- |
| [Adamson](../Tools/Omics/PerturbationPrediction/Adamson/EXPERIMENTAL_ROLES.md) | 50,440 selected cells and a fixed unseen-target protocol | Primary control/construct assignments and reconciliation of 94 observed guide groups versus the paper's 93-guide roster. |
| [Parse IFN-beta](../Tools/Omics/PerturbationPrediction/ParseIFNB/FEATURE_IDENTITY.md) | Complete counts, donor aggregates and DE | Intervention dose/reagent identity and an explicit feature contract: 409 duration-panel symbols are absent by exact name; naming candidates are not verified replacements. |
| [GSE181897](../Tools/Omics/PerturbationPrediction/GSE181897/RESULTS.md) | Primary author condition mapping; all 124 predictions, native reconstruction and independent scoring complete | HIRISA mean primary PASS, Kang mean primary FAIL; uncertainty remains deficient. New development on these now-scored outcomes needs a new validation cohort. |

Adamson and Parse remain unfitted and unscored for their proposed prediction tests.
Resolve their roles from experimental records, not RNA patterns. Freeze any
justified protocol revision before prediction and disclose earlier outcome
inspection, including Parse DE. Further storage improvements do not clear these
remaining admission gates.

The latest [native legacy H5AD check](../Tools/Omics/H5AD/Projection/README.md#original-legacy-kang-2026-09-11)
preserves every original Kang cell, gene, annotation and embedding, including
its category definitions. This makes the original information available without
the earlier selected-metadata conversion. It is an interoperability result;
no new prediction fit, held-out outcome or biological generalization is claimed.
The subsequent [complete analytical route](../Tools/Omics/H5AD/Projection/LEGACY_COUNT_ROUTE.md)
also verifies every native count record and donor/condition/cell-type aggregate,
with source RNA totals matching all cells. That closes the input-to-aggregation
handoff; the predictive conclusions and remaining requirements below still apply.

The [profile-guided count-writer qualification](../Tools/Omics/CountStore/Metal/Profile/README.md)
checks the same complete source information under a revised storage implementation.
Exact count/normalization bytes and bounded memory-safety checks protect that
input-processing boundary. They add no held-out biological observations: the
remaining decision is whether frozen predictors outperform simple baselines in
an independently identified experiment, with uncertainty and failures reported.

The [complete Replogle 2020 UPR prediction experiment](../Tools/Omics/PerturbationPrediction/Replogle2020/RESULTS.md)
now meets the fixed primary comparison in all five original technical gemgroups.
All 150 held-target predictions and native replays completed. Equal-target
all-gene RMSE is 1.41–2.12% below the training mean, but 34/150 folds remain worse
than no change and 37/150 worse than the mean. The [identity audit](../Tools/Omics/PerturbationPrediction/Replogle2020/IDENTITIES.md)
resolved nominal gene descriptors using a reproduced guide table and reference
sequence evidence; it does not prove full guide identity or off-target specificity.
This is separate data collection in the same broad K562/UPR setting, with shared
investigators and earlier target selection, not independent laboratory or tissue
validation. Technical gemgroups are not biological replicates.

The separate [GSE181897 external B-cell cohort](../Tools/Omics/PerturbationPrediction/GSE181897/README.md)
has passed complete source-count validation and native B-lineage aggregation:
136,142 source cells, 64 donor IDs, and all 34,287,682 selected RNA records are
checked. The source labels all features as gene expression, but its genome
field separates RNA from antibody counts and reproduces both author modality
totals exactly. Historical primary author notebooks now resolve B as IFN-beta
and C as control. The [complete frozen prediction test](../Tools/Omics/PerturbationPrediction/GSE181897/RESULTS.md)
uses all 62 eligible pairs, 4,938 source B-lineage cells and 11,800 exact shared
genes, retaining unmatched donors 5/23 and 84 absent panel genes. All 124 native
predictions reconstruct; all 496 estimate vectors and interval bounds match
independent calculations, and repeated scoring is byte-identical.

HIRISA-trained mean improves equal-donor RMSE by 5.62%, passing the frozen 5%
criterion and beating no-change in every donor. Kang-trained mean worsens RMSE
by 0.89% and loses in 41/62 donors. Kang ridge passes its secondary comparison;
HIRISA ridge loses to its mean in all 62 donors. Nominal 95% treated-expression
coverage averages 92.39% for Kang and 35.54% for HIRISA, with donor ranges
79.14–97.47% and 30.95–42.84%. The HIRISA intervals are much narrower and severely
under-cover. This adds an external RNA point-prediction success for one frozen
training origin, not calibrated uncertainty, universal transfer, verified
participant independence or a measured phenotype model.

The [retrospective cell-sampling diagnosis](../Tools/Omics/PerturbationPrediction/GSE181897/Uncertainty/README.md)
then resamples all 4,938 admitted original cells within the 124 donor/condition
strata. Observed sampling variability is substantial, while zero-count genes
remain poorly covered. Outcome-informed widening using both control and treated
cells raises HIRISA treated coverage to 94.25%; it is not available at query time
and is not a new predictor. The next uncertainty extension needs an explicit
count-based observation model, source sampling information and zero-count
support, followed by new validation. The original external interval failure stands.

The [complete GSE226572 external experiment](../Tools/Omics/PerturbationPrediction/GSE226572/README.md)
now evaluates native whole-population predictions across all three new donors
and all 18 released IFN-beta donor/time profiles. Every one of the 24 source
libraries was checked; frozen initial QC admits 126,633 cells. Training uses all
24,673 admitted Kang cells, with no query-response fitting. Mean response reduces
equal-donor, equal-within-donor-time RMSE by 2.00%, below the predeclared 5% target,
and is worse than no change in 8/18 cases. Ridge also fails; nominal 95% interval
coverage ranges from 85.37% to 95.45%. All native estimates, bounds and repeated
scores pass numerical checks. This separately collected, whole-population RNA
test is distinct from the B-cell endpoint and does not reproduce the author's
later curated cell selection. The unchanged six-hour response has no duration
covariate; differing time, health, culture, dose, assay, QC and composition limit
interpretation. Any repair using these now-scored outcomes needs new validation.

The subsequent [native exposure-duration model](../Tools/Omics/PerturbationPrediction/Duration/README.md)
now completes all three donor-held-out folds and all 18 outcomes in GSE226572.
Averaging training-donor curves interpolated in log1p(hours) reduces RMSE by
**55.71% versus no change** and **24.83% versus a matched time-invariant mean**,
passing its predeclared development gate. Duration mean beats no change in all
18 cases but loses to the matched mean twice; context ridge fails its secondary
gate. All native reconstructions, independent point/bound checks and repeated
scores pass. Nominal 95% treated-expression coverage is 93.59–96.17%, with broad
intervals from only two training donors per fold. This is retrospective model
development on inspected data, not independent validation or general calibration.
It identifies useful timing information within this study while preserving the
original external fixed-response failure and the missing phenotype boundary.

The next [Parse IFN-beta source admission](../Tools/Omics/PerturbationPrediction/ParseIFNB/README.md)
selects all 725,031 released IFN-beta/PBS cells across 12 source donor IDs
and retains every source RNA feature. All twelve native donor ingestions and
independent cell/aggregate checks now pass, covering 1,373,870,697 records and
3,070,817,047 retained counts. [All twelve native source replays and complete offline artifact restoration now pass](../Tools/Omics/PerturbationPrediction/ParseIFNB/COUNT_RESULTS.md). The first monolithic
attempt stopped on upstream HTTP 500 and published no complete result; bounded
retries and saved donor receipts protect the replacement run. No predictor has
been fitted or scored. The source historical
`gene_count` differs from retained X cardinality for 30,634 cells. Exact symbol
correspondence also leaves 409 of the duration model's response features absent,
and its dose/reagent identity is unresolved. A subsequent [HGNC nomenclature
audit](../Tools/Omics/PerturbationPrediction/ParseIFNB/FEATURE_IDENTITY.md) identifies
271 unique candidates but leaves 138 unresolved/ambiguous; it changes neither
the panel nor the prediction gate. Neither QC agreement nor model
compatibility is forced by excluding cells, filling absent genes with zero or
silently changing the panel. The count-stream owner checks arithmetic and
reconstruction; the Python/HDF5 adapter owns remote extraction and source ranges.
Its receipts do not claim a hash or count validation of the entire 227 GB file.

The [complete Parse file-backed differential-expression qualification](../Tools/Omics/CountStore/Expression/README.md)
now connects the count bundle to the existing paired negative-binomial model.
All 725,031 cells and 24 donor/condition aggregates are included. The new owner
exactly reproduces the prior native statistical report apart from explicit
membership encoding; all six edgeR, limma-voom and DESeq2 comparisons complete.
There are 33,899 tested features, 5,946 rank-deficient-support exclusions and 507
low-count exclusions. Effect concordance with external methods is a useful
statistical check, but does not establish false-discovery calibration, a causal
cell-intrinsic effect or an unseen-outcome prediction. The [incremental report writer](../Tools/Omics/CountStore/Expression/Memory/README.md)
now lowers the same complete analysis from 970.4 to 202.3 MiB native peak RSS
while preserving every report byte. Per-feature model state remains resident;
this storage change adds no biological outcome evidence. No Parse duration predictor has been fitted or scored;
the panel/dose gates above and the prior transfer failures remain in force.

The [native Visium import](../Tools/Omics/Multimodal/Visium/README.md) now makes
original spatial RNA counts and pixel positions available without an H5MU
preparation step. Every count and coordinate in the complete 4,039-spot source
matches independent readers and the earlier qualified interchange product.
These measured spatial inputs do not themselves supply a validated spatial
response model, deconvolution, tissue function or variant-to-phenotype endpoint.

The [native balanced reference classifier](../Tools/Omics/ReferenceMapping/Logistic/README.md)
now predicts candidate source labels for all 8,569 Baron query cells, reproducing
the previous external classifier. Macro-F1 improves over kNN in every donor,
with lower overall accuracy in two and retained rare-class misses. This advances
the annotation step; calibrated confidence, authoritative biological identities
and independent context transfer remain unqualified. It adds no perturbation
response, tissue-function or clinical outcome evidence.

The [complete Kang–Ding annotation-transfer experiment](../Tools/Omics/ReferenceMapping/CrossStudy/README.md)
now retains all 68,704 query cells and each original RNA-library denominator using
an explicit shared gene panel. Numerical reconstruction and independent checks
pass, but both declared coarse-family targets fail: macro-F1 is 0.613 and 0.666,
and every source-labelled megakaryocyte is missed in each direction (800 and 132).
Ding labels cover only 66.80% of its query cells; original fine taxonomies differ.
This is evidence against generalizing the earlier within-study annotation success,
not evidence of reliable biological identity or outcome prediction across studies.

The [source-label and frozen-model diagnosis](../Tools/Omics/ReferenceMapping/CrossStudy/Diagnosis/README.md)
checks all 24,673 original Kang labels and all 68,704 query cells. The preparation
preserved the original labels exactly. Kang's literal rare class has strong
HBB/HBA1/HBA2 expression and little PF4/PPBP/ITGA2B/GP9/TUBB1 expression; Ding's
same-named class shows the opposite pattern. All diagnostic genes were retained
and both classifiers recognize their own training class almost perfectly. Exact
gene-space reconstruction traces the failed margins without refitting. These
post-result observations raise a source-label-validity concern; they do not
establish replacement identities, revised success or a complete causal account.
Both original failures remain, and cross-study label equivalence needs independent
support before a new annotation-transfer qualification.

## What information is sufficient for the implemented predictors?

| Question | Information available before prediction | Output and present evidence |
| --- | --- | --- |
| How will a new donor respond to a known treatment? | Paired control/treated count profiles from training donors; the new donor's control counts; matching treatment, organism, population and feature identities. | Native no-change, mean, median and context-ridge estimates of population-average treated log1p(CPM). Complete Kang, Hagai and HIRISA held-out donor experiments show conditional predictive signal, with failures retained. |
| How does a known treatment's RNA response depend on exposure duration? | Matched training-donor control and at least two timed treatment profiles, explicit exposure hours and matching intervention/population identities, plus the held-out donor's control RNA. | Native donor-curve interpolation within the observed time range. Complete three-donor GSE226572 development achieves 24.83% lower RMSE than matched time-invariant mean; independent temporal transfer and extrapolation remain unqualified. |
| What happens when two observed targets are perturbed together? | Control and single-target counts for both constituents in the same experimental context; the pair's target identities. | Six native composition baselines evaluated on all 131 held-out Norman pairs. This tests expression composition; it does not identify genetic interactions. |
| What happens when an unobserved target is perturbed? | Control and other single-target counts in the same context; independently supplied target identities and GO terms for training and query targets. | Native target-kernel prototype: 105 Norman development folds (101 supported), followed by 150 Replogle folds (30 targets × 5 technical groups). Replogle meets the fixed primary criterion in all five groups with modest gains and target-level failures; broad unseen-target transfer remains unqualified. |
| Can a known response transfer between enriched preparations and PBMCs? | Other-donor paired counts from the training preparation; query-preparation control counts; exact author lineage annotations, donor and feature identities. | All 120 HIRISA folds completed. Cross ridge beats both no-change and cross mean in 3/12 contrasts; this is an inspected-study, annotation-conditioned endpoint with preparation/batch confounding. |
| Can the IFN-beta response transfer between Kang and HIRISA? | Paired counts from the other study, query-donor control counts and an explicit shared-gene mapping, preserving each library's full measured-feature denominator. | All 26 cross/within folds completed over 11,884 matched genes. Cross ridge fails both directional primary comparisons and is worse than cross mean for all 13 donors. This is a combined health, duration, preparation and assay shift on reused studies. |
| What will happen in a new tissue, species, disease state or patient? | Would require a validated transfer model and measured outcomes in that destination context. | No qualifying result in this evidence set. Existing point estimates cannot be promoted to those outcomes. |
| Can a DNA variant predict a cellular or tissue phenotype? | Would require assembly/allele-resolved regulatory evidence and validated connections through RNA, proteins, mechanisms and phenotype. | [AlphaGenome integration work](AlphaGenomeAtlas.md) provides a direction and separate evidence interface; a validated end-to-end phenotype predictor is not established. |

The expression outputs are compositional pseudobulk point estimates, not absolute
molecule counts or individual-cell response distributions. Negative predicted
log-expression is clipped to zero; implied CPM totals are retained rather than
silently renormalized. Optional normal-model mean-response intervals are now
implemented and assessed below; they are not generally calibrated. These models
still supply no experimentally validated RNA-to-phenotype mapping.

## Decision before using a prediction

### Available information does not imply a validated outcome

There are three distinct questions: whether the original observations and their
identities are available, whether an implemented model accepts those inputs for
the requested endpoint, and whether its predictions have passed an appropriate
held-out comparison. A positive answer to the first question cannot substitute
for the other two. For example, the complete Parse IFN-beta ingestion supplies
725,031 cells, but unresolved intervention and feature identities still prevent
admitting that cohort to the unchanged duration-model protocol. No Parse outcome
has been fitted or scored.

The complete donor replay and artifact restoration have now passed. The next
prediction step is to resolve the intervention and feature contract before
freezing an independent prediction comparison. Any explicitly revised panel or
model requires a new protocol and an account of which outcomes were already
inspected; the 271 HGNC naming candidates do not silently repair the original
409 absent symbols. Report every donor against no-change and the matched mean,
alongside uncertainty coverage, interval width and failures. A failed comparison
remains a result rather than a reason to change the acceptance threshold.

The [file-backed cell-axis and count consumer](../Tools/Omics/CountStore/CellAxis/README.md)
now address resident cell identities, QC and membership arrays for count aggregation.
Fourteen native tests pass, and all 725,031 original Parse cell identities and
declared matrix totals pass import and reopen at 62.6 MiB native peak RSS.
The [complete paired count ingestion, source replay and memory comparison](../Tools/Omics/CountStore/CellAxis/COUNT_RESULTS.md)
now pass: file-backed ingestion/replay peak at 165.6/181.8 MiB versus
1,005.0/1,024.7 MiB for the resident owner. These completed count-path results do
not establish that graphs or prediction models execute out of core, nor do they
add a biological prediction result. Atlas supplies a separate
variant-to-molecular hypothesis route; it supplies neither missing Parse treatment
metadata nor experimental RNA-to-phenotype outcomes.

### Applying the present evidence

The available evidence supports a **conditional research estimate of average RNA
response** when the experiment matches an implemented predictor's input contract.
It does not establish that any particular new query will beat a simple baseline.
The aggregate Replogle improvement cannot identify the failing targets in advance;
the HIRISA ridge results likewise do not justify replacing the mean response
across all populations and treatments.

For a proposed use, first specify the observable, units, intervention, population,
and time of measurement. Supply the control and training observations required by
the input table. Keep no-change and training-mean estimates alongside the proposed
model, its provenance and unsupported inputs. If the requested endpoint is tissue
function, survival or treatment benefit, the current RNA output leaves that
endpoint unanswered. Additional input annotations alone do not close the missing
outcome model or validation.

### Uncertainty has a separate owner and evidence requirement

The native [single-cell response model](../Sources/NumiVivoKit/Omics/VivoPerturbation.swift)
stores mean/median responses and ridge coefficients. Its optional
`donorResponseIntervalCoverage` now adds a separate normal-model future-donor
interval around the mean response; it does not supply context-ridge uncertainty.
The [empirical assessment](../Tools/Omics/PerturbationPrediction/Intervals/README.md)
retains coverage, width, unavailable genes and the failed context-transfer
assumption rather than promoting nominal coverage. The existing
[Bayesian target-engagement predictor](../Sources/NumiVivoKit/Calibration/VivoTargetPosteriorPrediction.swift)
instead propagates kinetic posterior particles into occupancy observables at
specified times, with a declared assay model. Its pointwise intervals, noise
assumptions and failure handling belong to that model. They cannot be attached
to gene-expression estimates without an explicit statistical model connecting
those quantities and separate empirical validation.

### Next prediction experiment

The next milestone is an independent biological-context test of a frozen
expression predictor, with uncertainty evaluated on held-out biological
replicates. Before inspecting outcome scores:

1. Admit a cohort using primary intervention, control, donor, context and feature
   identities. Define which observations may be available at query time; do not
   require treated-query annotations that a prospective user would lack.
2. Freeze training, any calibration data, and final test donors without donor
   overlap. Record prior inspection and study/model-development overlap. Keep
   technical groups nested under their biological replicate rather than counting
   cells, genes or overlapping folds as independent donors.
3. Fix the endpoint, full-feature denominator, model, simple baselines, primary
   comparison and minimum useful improvement. Preserve unsupported cases and
   failures. Report per-donor and per-context results as well as the aggregate.
4. If adding intervals, declare their statistical assumptions, nominal coverage,
   whether coverage is per gene or simultaneous, and how calibration uses only
   non-test donors. Assess held-out coverage together with interval width and
   missingness; a nominal probability or a wide interval alone does not establish
   useful uncertainty. Report insufficient replication explicitly.
5. Freeze native predictions before scoring, retain source/executable identities,
   verify numerical reconstruction independently, then report the biological
   result even when it fails. Development after test inspection requires a new
   validation cohort before making a stronger generalization claim.

These are acceptance requirements for stronger generalization claims. Completed
experiments below retain their original protocols and limits. The new nominal
interval implementation does not close the independent-calibration requirement.

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

The [retrospective donor-score sensitivity check](../Tools/Omics/Benchmarks/HIRISA/DONOR_SCORE_SENSITIVITY.md)
reconstructs all published scores and omits each scoring donor in turn, without
refitting. Thirteen of the fourteen ridge wins over no-change survive every
omission. NK IFN-L1 changes sign; the already negative Monocyte IFN-L1 comparison
is also sensitive. Only three of four ridge wins over the training mean survive
every omission: Monocyte IFNa, IFNb and IFNg. NK IFNa is sensitive. These are
aggregation diagnostics on reused evidence, not confidence intervals or new
independent prediction experiments.

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

### Known treatment across preparations

The [complete HIRISA transfer experiment](../Tools/Omics/Benchmarks/HIRISA/CONTEXT_TRANSFER.md)
withholds each query donor from both preparations and retains six separate author
lineages. All **120 frozen folds completed**: sixty cross-preparation predictions and
sixty matched within-preparation references, spanning 705,365 selected cells and
all 18,082 source genes. All six native count bundles and prediction batches
pass independent numerical checks and native replay. Every output was frozen
before scoring, and two scoring executions produced identical result bytes.

All three learned methods beat no-change in **all twelve cross-preparation
contrast means**. Ridge meets the fixed
primary gate—lower all-gene RMSE than both no-change and cross-preparation mean—in
**3/12 contrasts**. Cross-preparation ridge is worse than matched within-preparation
ridge in all twelve contrasts, showing a consistent transfer penalty. No model
was selected or tuned from these results.

The same donor-score audit retains all 120 cross/within folds. All twelve
cross-ridge improvements over no-change and all twelve cross-versus-within ridge
penalties survive every single-donor score omission. Nevertheless, cross ridge
is worse than no-change in four individual donor folds and beats cross mean in
only 21/60 individual folds. All three contrast-mean wins over cross mean survive
omission. Shared training donors and reused contexts prevent treating these
folds as independent biological experiments.

The protocol preceded these fits and PBMC scores, but followed other inspected
HIRISA results. Treated-cell annotations define outcome strata, and preparation,
batch and culture composition vary together. This supports measured conditional
average RNA prediction across these combined contexts; it does not establish
independent-study replication, prospective cell identity or a causal preparation
effect. All twelve direction/lineage outcomes and every fold remain reported.

### Known treatment across studies: failed transfer comparison

The [Kang–HIRISA experiment](../Tools/Omics/PerturbationPrediction/CrossStudyIFNB/README.md)
uses all eight Kang and five HIRISA B-cell donor pairs, with 11,884 exact, unique
source-symbol matches. Each library retains all its own measured genes in the
normalization denominator. Native fits/replays and independent reconstruction
pass for all 26 cross/within folds and all 104 prediction vectors; repeated
scoring is identical. Source-qualified pseudobulk transport is explicitly
distinguished from re-running raw single-cell ingestion.

| Training → query | No change | Cross mean | Cross ridge | Within ridge |
| --- | ---: | ---: | ---: | ---: |
| HIRISA → Kang | 1.220670 | 1.135572 | 1.178674 | 1.132625 |
| Kang → HIRISA | 0.341591 | 0.521068 | 0.596829 | 0.107549 |

These are equal-donor mean response RMSEs over the complete shared panel. Ridge
must beat both no-change and cross mean to pass the fixed primary comparison:
**neither direction passes**. Cross ridge is worse than cross mean for every
donor, and worse than no-change for all five HIRISA query donors. HIRISA's simple
mean response has some transfer signal toward Kang, but no predictor is promoted.

The query cohorts differ in disease status, stimulation duration, preparation,
cell annotation and measurement chemistry. These results do not isolate the
cause of the transfer failure. Both studies were previously inspected; this is
a new cross-study test on reused observations, not untouched external validation.
Matching gene symbols also does not prove a common reference annotation release
or assay equivalence. The negative result strengthens the present limit on
generalizing the donor-context predictor.

### Mean-response intervals: measured coverage is context dependent

The [native interval experiment](../Tools/Omics/PerturbationPrediction/Intervals/README.md)
adds nominal 95% future-donor intervals to the same 26 frozen Kang–HIRISA folds.
Training sample variances and Student-t quantiles determine bounds around the
mean response. Every fit/prediction passes replay; independent NumPy/SciPy checks
verify all 104 bound arrays, and the 104 existing point vectors remain exact.
Constant-response genes retain unavailable intervals rather than zero-width
certainty. Clipped and unclipped coverage are reported separately.

| Query study / mode | Mean treated coverage | Mean treated width | Gene availability |
| --- | ---: | ---: | ---: |
| Kang / within | 92.94% | 4.230723 | 87.65% |
| Kang / cross from HIRISA | 42.39% | 0.279829 | 99.98% |
| HIRISA / within | 94.70% | 0.467843 | 99.97% |
| HIRISA / cross from Kang | 99.21% | 5.152851 | 88.57% |

Values are equal-donor averages among available genes; widths use log1p(CPM).
Within-Kang donor coverage ranges from 81.45% to 98.92%, so the mean conceals
substantial variation. Cross-to-Kang raw-response coverage is only 33.08%; clipping
raises it to 42.39% through the zero boundary. Reverse transfer covers broadly
with wide intervals. Neither high aggregate coverage nor nominal 95% establishes
useful calibration. Normal exchangeable responses, a fixed query control and
pointwise coverage remain assumptions; these reused small-donor cohorts do not
qualify simultaneous gene coverage or general uncertainty.

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

### Unseen targets: development and separately collected validation

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

The subsequent [Replogle validation](../Tools/Omics/PerturbationPrediction/Replogle2020/RESULTS.md)
uses the unchanged lambda-one native model and all source RNA genes. All five
gemgroups meet the predeclared all-gene criterion against the all-training mean,
supported-training mean and one-position shuffled response. Both means coincide
because all 30 targets have usable direct GO descriptors. Gains over the mean
range from 1.41% to 2.12%; no-change regressions affect 34/150 folds and mean
regressions affect 37/150. All 750 vectors pass independent checks, and repeated
scoring is identical. No model or favorable group was selected after scores.

This strengthens the evidence beyond Norman's reused development data. It tests
an independently collected experiment with shared investigators, K562 systems
and prior UPR target selection. It does not transfer Norman-trained coefficients,
provide temporal independence of current GO knowledge, or establish prospective
biological contexts. Nominal target identity relies partly on FBA's reproduced
sequence table; the original full guide supplement was not retrieved.

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
controls remain insufficient, and complete preservation remains unqualified.

The subsequent [annotation-retention diagnostic](../Tools/Omics/Benchmarks/HIRISA/ANNOTATION_RETENTION.md)
evaluates every original cell across 114 donor folds, 23 preparation/treatment
strata and 31 exact author labels. Native integration exceeds recall-loss margins
in 37/146 supported, control-sensitive label/stratum comparisons, including
10/29 rare comparisons. Harmony seeds 7/19/41 fail 41/39/38 comparisons. All five
representations match independent per-cell SVD fits and exact query confusion.
Another 308 comparisons lack sufficient support and seventeen more lack control
sensitivity. These are annotation recoverability failures under a fixed decoder;
author predictions are not biological truth. No candidate qualifies complete
preservation, and native multi-seed robustness remains open. Full-cohort [clustering publication, native replay and independent
partition checks](../Tools/Omics/Benchmarks/HIRISA/FULL_CLUSTERING_RESULTS.md) now
pass for all original cells; graph communities are not validated cell types.

A subsequent [protected-stratum regression experiment](../Tools/Omics/Benchmarks/HIRISA/PROTECTED_REGRESSION.md)
conditions donor effects on known preparation/treatment strata with the same
frozen native memberships. It retains 37/146 annotation failures, resolving four
and introducing four; rare failures fall from 10 to 8. Response and program
gates are unchanged. The modest average recall improvement does not repair
complete preservation, and no native solver change is promoted.

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
   decoder/Kang NK-recall failures. Address the measured HIRISA annotation/rare-label
   retention failures, retain all insufficient strata, and complete neighborhood
   and independent biological validation.
   The [compact native program artifact](../Tools/Omics/Benchmarks/HIRISA/NATIVE_PROGRAM_RESULTS.md)
   now passes full-cohort publication, replay and independent comparison.
3. Resolve Adamson controls and its 94-versus-93 guide roster from primary records, then execute the
   frozen independent-study target-prediction protocol with coverage, all
   failures and matched simple/shuffled baselines. Do not tune it on test scores.
   The complete Replogle fixed-model experiment now passes its declared primary
   comparison in all five technical groups. Preserve all 150 folds and failures,
   and extend validation to an independently selected target panel and laboratory
   or biological context rather than retuning on these results. Adamson's gate
   remains unchanged.
4. Establish transfer with prospective strata and measured outcomes in an
   independent study. The [complete GSE226572 test](../Tools/Omics/PerturbationPrediction/GSE226572/README.md)
   now adds a separately collected experiment with all donor/time outcomes;
   it fails its 5% mean-improvement and ridge targets. The subsequent
   [native duration model](../Tools/Omics/PerturbationPrediction/Duration/README.md)
   passes the complete donor-held-out development gate, while its context ridge
   does not improve on duration mean. Freeze this model for a new independent
   temporal experiment with matching dose, timing and population metadata; preserve
   the fixed-model failure and report every held-out outcome. The new [Kang–HIRISA comparison](../Tools/Omics/PerturbationPrediction/CrossStudyIFNB/README.md)
   fails both directional primary gates on reused studies; it does not close this
   requirement. The [completed HIRISA preparation-transfer experiment](../Tools/Omics/Benchmarks/HIRISA/CONTEXT_TRANSFER.md)'s
   sixty cross-preparation folds and sixty matched references now pass replay
   and numerical checks; limited donor replication and context confounding remain.
5. Establish predictive interval coverage and useful improvements over simple
   baselines on independent biological replicates before promoting a predictor.
   The new normal-model mean-response interval is implemented, but its
   [reused-cohort assessment](../Tools/Omics/PerturbationPrediction/Intervals/README.md)
   exposes undercoverage, broad intervals and missingness; this gate remains open.
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
