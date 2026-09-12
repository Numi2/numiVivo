# NumiVivo

**From atoms to biological behavior. Built for Apple silicon.**

NumiVivo is a research platform for molecular simulation and programmable biology. It brings molecular structures, GPU molecular dynamics, native electronic structure, reaction networks, single-cell analysis, expression-response prediction, and reproducible experiments into one Swift, C++ and Metal codebase.

The ambition is to follow a molecular change across scales: how a structure moves, where electrons interact, which reactions become possible, and how those reactions affect a cell or physiological model. Every transition should preserve the quantities, assumptions and evidence needed to understand the result.

[Get started](#get-started) · [Explore the capabilities](#explore-the-capabilities) · [Examples](#choose-an-experiment) · [Documentation](Documentation/README.md) · [Implementation status](Documentation/CAPABILITIES.md) · [Completion roadmap](Documentation/COMPLETION_ROADMAP.md)

> **Research software, under active development.** The capabilities below describe source implementations and their intended workflows—not a fully qualified release. Apple package integration, GPU numerical behavior and performance require qualification. The [capability map](Documentation/CAPABILITIES.md) separates implemented methods, current restrictions and planned work.

The [source-bound atlas export](Tools/Omics/H5AD/Annotation/Atlas/Clustering/README.md) now carries all 1,612,594 native graph-community labels into a backed AnnData-readable H5AD, with complete cell-identity and metadata checks. This is interoperability evidence, not biological annotation validation.

The [integration design guard](Tools/Omics/Reduction/INTEGRATION.md#optional-protected-sample-groups) can reject correction designs confounded with explicitly protected sample groups. Admission does not establish preservation of biological signals.

The [full-HIRISA rigid-correction experiment](Tools/Omics/Benchmarks/HIRISA/RigidProjection/README.md) preserves within-donor distances but still fails 31/146 sensitive annotation comparisons. No new integration method is promoted.

The [translation-only ablation](Tools/Omics/Benchmarks/HIRISA/TranslationOnly/README.md) reduces annotation failures to 22/146 but fails a B-cell IFNa response comparison and leaves more donor scatter. Fewer failures do not establish adequate integration; no method is promoted.

The [matched-control donor-shift experiment](Tools/Omics/Benchmarks/HIRISA/MatchedControlShift/README.md) estimates shifts from 403,449 control cells but still fails 34/146 annotation comparisons and one response contrast. Control-only shift estimation does not resolve the preservation gap.



The [exact CPU matcher is now integrated into the Swift MNN owner](Tools/Omics/Reduction/ExactNativeMNN/README.md). Actual execution reproduces scores, anchors and reports byte-for-byte across all 82,567 Hagai/Kang/Ding cells; sanitizer and cancellation checks pass. The earlier [biological benchmark](Tools/Omics/Reduction/ExactCPUMNN/README.md) passes every evaluable gate, with four Kang classifier strata still unavailable. Million-cell and independent biological qualification remain open.

## Biological prediction: current evidence

A [training-only response-shrinkage candidate](Tools/Omics/PerturbationPrediction/StudyContextKernel/ResponseShrinkage/README.md) improves HIRISA but still fails Kang and GSE181897 transfer gates. All 75 donors pass numerical verification; the candidate remains unpromoted.

The [completed Parse context prediction test](Tools/Omics/PerturbationPrediction/ParseIFNB/ContextEvaluation/README.md) fails: mean donor RMSE is 20.93% worse than training-mean response and 32.14% worse than no change. All twelve donors lose to training mean. All 139,200 native predictions and every donor score pass independent numerical checks. This supersedes earlier statements that Parse prediction is pending; the candidate remains unpromoted, and participant independence and matched exposure remain unverified.


The [paired RNA/ATAC product workflow](Tools/Omics/Multimodal/PairedWorkflow/README.md) now runs original 10x counts through native RNA PCA, ATAC TF-IDF/LSI, cell-specific modality weights and weighted graphs. All stages match the separately qualified 2,711-nucleus benchmark, and complete source reconstruction passes. Execution remains bounded and resident. The [paired reconstruction diagnostics](Tools/Omics/Multimodal/ATACPreservation/README.md) fail both modalities’ combined no-worsening rules: primary RNA/ATAC reconstruction worsens versus modality-only baselines, despite small fuzzy-graph improvements. Biological preservation, unseen-context prediction and regulatory validity remain unproven. A [variance-retaining RNA candidate](Tools/Omics/Multimodal/RNAVarianceProfile/README.md) also fails both preservation gates against fixed and matched baselines; the product default is unchanged.

The [exact Parse-compatible panel evaluation](Tools/Omics/PerturbationPrediction/StudyContextKernel/ExactParsePanel/README.md) reran all 75 development donors with a frozen 11,600-feature contract. HIRISA passes the existing transfer gate; Kang and GSE181897 fail. That development experiment did not score Parse outcomes; the subsequent completed Parse test above fails, and the candidate remains unpromoted.

The [CellTypist CLI](Tools/Omics/ReferenceMapping/CellTypist/README.md) reproduces imported-model inference on 267 real Parse cells (529,404 count records), from canonical streams or native AnnData CSR/CSC files. H5AD runs retain the complete source and exact cell/feature identities; both sparse layouts produce identical output. Invalid inputs leave no completed bundle. This scoped verification does not establish held-out annotation accuracy or general biological prediction.

**Available data supports conditional estimates of average RNA responses;
reliable prediction of general biological outcomes is not established.**
NumiVivo can generate these estimates when the required training measurements,
control profile and biological identities are available. Held-out experiments
show useful results in some settings and failures in others; having compatible
inputs does not establish accuracy for a new query.

The [native multigene context model](Tools/Omics/PerturbationPrediction/StudyContextKernel/README.md) now uses the full untreated expression profile. It improves held-out HIRISA RMSE by 21.1% over the training-mean response, but is slightly worse than that baseline on Kang and GSE181897. All 75 donors and 885,000 predictions are verified. Reliable across-study transfer remains unestablished; these are reused development cohorts, not a fresh external test.

| Requested outcome | Decision from the available evidence |
| --- | --- |
| Average RNA response with matched training and untreated query measurements | **Supported as a conditional research estimate.** Performance depends on the treatment and population; retain no-change and training-mean baselines. |
| RNA response to an unseen gene target with training perturbations, controls and gene annotations | **Limited validation.** Fixed GO prediction improves mean RMSE by 1.41–2.12% across five Replogle technical groups, but loses to the mean in 37/150 target/group folds; these are not independent biological replications. |
| Exposure-dependent RNA response | **Development result.** Duration mean improves RMSE by 24.83% over the matched time-invariant mean in three donor holdouts; independent validation remains open. |
| RNA response in another study | **Unreliable transfer.** HIRISA-trained mean passes the GSE181897 target, but Kang-trained mean fails. The completed Parse context test is 20.93% worse than training mean and 32.14% worse than no change; reliable uncertainty remains unestablished. |
| Molecular effects of a DNA substitution | **External hypothesis source.** AlphaGenome Atlas can contribute variant evidence; its results do not qualify a NumiVivo phenotype prediction. |
| Individual-cell behavior, tissue function, disease progression or treatment benefit | **Not established.** These endpoints need their own outcome models and experimental validation. |

Successful count ingestion, numerical replay and large-cohort processing establish
that information can be processed correctly. They do not establish that an
unmeasured biological outcome can be predicted accurately. The
[input-to-outcome decision](Documentation/BiologicalPrediction.md#available-information-does-not-imply-a-validated-outcome)
separates these requirements and identifies the next work.

For a concrete prediction, specify the measured endpoint and units, population,
intervention and exposure time, then identify which training and control
measurements would be available before the outcome. The
[current completion and admission gates](Documentation/BiologicalPrediction.md#current-completion-and-admission-gates)
separate finished execution work from the missing biological evidence.

In 79 HIRISA held-out donor folds, context ridge beats
no-change in 14 of 16 contrasts, but beats the simpler training-mean response in
only four. Two contrasts are worse than no-change under every learned baseline.
A [retrospective donor-score check](Tools/Omics/Benchmarks/HIRISA/DONOR_SCORE_SENSITIVITY.md)
finds that 13 of those 14 no-change wins, and three of the four mean-baseline
wins, survive every single-donor score omission. This does not establish
predictive uncertainty or independent-study validation.
Native combination models also predict held-out pairs of observed targets;
the unseen-target GO prototype improves the mean baseline by 0.58% on reused
Norman development data. A separately collected [Replogle UPR validation](Tools/Omics/PerturbationPrediction/Replogle2020/RESULTS.md)
now meets the fixed primary criterion in all five technical gemgroups, with
1.41–2.12% lower all-gene RMSE than the training mean. It remains worse than
no-change in 34/150 folds; this does not establish general target or tissue transfer.

The [biological prediction assessment](Documentation/BiologicalPrediction.md)
explains the required inputs, complete positive and negative results, native
commands and remaining validation. These are expression point estimates;
reliable prediction of new tissues, disease outcomes, individual-cell responses
or variant-to-phenotype effects is not established.

The next prediction milestone is useful, reliable improvement in independent
biological contexts, with endpoints, donor splits, simple baselines and uncertainty
assessment fixed before scoring. [Acceptance requirements](Documentation/BiologicalPrediction.md#decision-before-using-a-prediction)
explain what would justify extending the present claims. Existing Bayesian
kinetic-model intervals do not provide uncertainty for the single-cell predictors.

The independent [Adamson validation](Tools/Omics/PerturbationPrediction/Adamson/EXPERIMENTAL_ROLES.md)
has verified 50,440 selected cells, but still needs primary control assignments
and reconciliation of its guide roster before fitting. [AlphaGenome Atlas](Documentation/AlphaGenomeAtlas.md)
can supply variant-level molecular hypotheses; its predictions do not establish
NumiVivo's downstream cellular or tissue outcomes.

A [Kang–HIRISA cross-study test](Tools/Omics/PerturbationPrediction/CrossStudyIFNB/README.md)
now completes all 26 cross/within folds over 11,884 shared genes. Context ridge
fails the primary comparison in both transfer directions and is worse than the
cross-study training mean for all thirteen donors. The studies differ in health
status, stimulation time, preparation and assay; their gene alignment does not
establish reliable biological-context transfer.

The [complete external PBMC duration-transfer experiment](Tools/Omics/PerturbationPrediction/GSE226572/README.md)
now checks all 24 GSE226572 libraries and 126,633 cells admitted by frozen initial
QC. Native predictions cover all three query donors and 18 released treatment-time
profiles. Mean response improves average RMSE by **2.00%**, missing the declared
5% target, and is worse than no change in 8/18 cases; ridge also fails. Nominal
95% interval coverage spans 85.37–95.45%. Numerical checks and repeated scoring
pass. A fixed six-hour response does not provide a validated temporal or general
biological-outcome model.

The subsequent [native duration-model development test](Tools/Omics/PerturbationPrediction/Duration/README.md)
uses all three donor holdouts and all 18 outcomes in that already inspected study.
Duration mean reduces RMSE by **55.71% versus no change** and **24.83% versus a
matched time-invariant mean**, passing the declared development gate. It loses to
the matched mean in 2/18 cases; context ridge fails to improve on duration mean.
All native reconstructions and independent checks pass. This supports conditional
RNA prediction when exposure time and matched training observations are available;
new independent validation is still required, and the earlier external failure stands.

The [Parse B-cell metadata audit](Tools/Omics/PerturbationPrediction/ParseIFNB/BCellAdmission/README.md) identifies 72,446 source-labeled B cells across all twelve paired donors. The original context model has 200 absent exact-name features. A versioned 11,600-feature panel now resolves exact source compatibility, but its renewed development evaluation fails transfer gates on two of three studies; the [primary dose/reagent is now known](Tools/Omics/PerturbationPrediction/ParseIFNB/BCellCounts/NEXT_STEPS.md#primary-dose-and-reagent-resolved), while cross-study exposure comparability and participant overlap remain unverified. The [completed Parse context prediction test](Tools/Omics/PerturbationPrediction/ParseIFNB/ContextEvaluation/README.md) now fails on this full frozen B-cell cohort: mean donor RMSE is 20.93% worse than training mean and 32.14% worse than no change. The [native B-cell count workflow](Tools/Omics/PerturbationPrediction/ParseIFNB/BCellCounts/README.md) now has all twelve native ingestions and source replays independently verified: 124,909,573 records and 253,128,870 total counts. The complete archive and its 13,990 members passed verification. This is count-handling evidence, not predictive validation.

The [Parse IFN-beta source cohort](Tools/Omics/PerturbationPrediction/ParseIFNB/README.md)
contains all 725,031 released IFN-beta/PBS cells from 12 donors. Native count
admission now passes all 1,373,870,697 selected records through resumable donor
streams; [all twelve native source replays and offline artifact restoration now pass](Tools/Omics/PerturbationPrediction/ParseIFNB/COUNT_RESULTS.md). This avoids a 227 GB local source copy. An initial upstream HTTP 500 failure and a
subsequently repaired temporary-memory defect are retained explicitly. Historical QC differs from retained matrix counts;
409 duration-panel symbols are absent by exact name. The source dose and reagent
are resolved; cross-study exposure equivalence remains unverified. These count
checks are separate from the failed B-cell context prediction test above.

The new [native file-backed cell axis and count consumer](Tools/Omics/CountStore/CellAxis/README.md)
remove resident cell identities, QC and membership arrays from the count path.
All 14 native tests pass, and importing and reopening the complete 725,031-cell
axis preserves every identity byte and declared matrix total at **62.6 MiB native
peak RSS**. [Complete paired ingestion and original-source replay now pass](Tools/Omics/CountStore/CellAxis/COUNT_RESULTS.md)
for all 1.37 billion records. File-backed ingestion/replay peak at **165.6/181.8 MiB**,
versus **1,005.0/1,024.7 MiB** for the resident owner on the same source.
This count-path memory improvement does not change the biological verdict.

The [file-backed count-to-DE analysis](Tools/Omics/CountStore/Expression/README.md)
now completes the full Parse cohort and exactly matches the existing native
statistical results. All six edgeR/limma/DESeq2 comparisons complete; shared-offset
effect correlations are 0.994, 0.860 and 0.990 respectively. Of 40,352 features,
33,899 are tested, 5,946 lack estimable support and 507 fail the count filter.
The [incremental report writer](Tools/Omics/CountStore/Expression/Memory/README.md)
reduces full-cohort analysis peak RSS from 970.4 to **202.3 MiB**, preserving
every report byte. These are observed RNA differences,
not a new held-out prediction or validation of tissue/clinical outcomes.

The [native interval assessment](Tools/Omics/PerturbationPrediction/Intervals/README.md)
now tests optional nominal 95% mean-response intervals on those same 26 folds.
Average treated-expression coverage is 92.94% / 94.70% within Kang / HIRISA,
versus 42.39% / 99.21% across studies. Width, missing intervals and donor-level
variation prevent interpreting these averages as general predictive calibration.

The [complete external GSE181897 B-cell prediction test](Tools/Omics/PerturbationPrediction/GSE181897/RESULTS.md)
now covers all 62 eligible donor pairs and 11,800 shared genes. Historical author
code resolves the treatment labels before fitting. HIRISA-trained mean reduces
average RMSE by **5.62%**, passing the frozen 5% target and beating no-change in
all 62 donors. Kang-trained mean worsens error by 0.89% and fails. All 124 native
predictions reconstruct; 496 estimate vectors match independent calculations.
HIRISA-trained nominal 95% treated-expression intervals cover only **35.54%** of
available features on average. This supports a bounded RNA point-prediction result,
with failed uncertainty and model-dependent transfer retained explicitly.
The [cell-sampling diagnosis](Tools/Omics/PerturbationPrediction/GSE181897/Uncertainty/README.md)
points to sampling and zero-count uncertainty; its outcome-informed widening is
a retrospective diagnostic, not a new prediction model.
The native [count observation posterior](Tools/Omics/CountObservation/README.md)
now models original cell counts and RNA depths, with nonzero uncertainty after
zero observations and separate variance for planned cell sampling. It is
conditional on its prior and dispersion. [Training calibration](Tools/Omics/CountObservation/Calibration/README.md)
now fits these from all 122,164 original Kang/HIRISA training cells and verifies
the path to available query-control posteriors. Unsupported genes remain explicit;
the joint response model and independent biological calibration have separate gates.
The [paired-donor diagnosis](Tools/Omics/CountObservation/Paired/README.md) now
checks every training gene and every donor omission. Separate noise corrections
produce invalid joint covariance for 7,489 Kang and 6,373 HIRISA genes; per-cell
and pseudobulk endpoint differences also matter. The new
[joint count-response model](Tools/Omics/CountObservation/Joint/README.md) fits
nonnegative paired rate distributions from all 122,164 training cells for the
previously selected 16-gene development panel. All 76 available finite-grid fits
and 4,712 conditional predictions pass independent numerical checks. However,
**17/19 available gene models fail grid refinement**; 13/32 origin/gene cases
remain unavailable. This implementation does not yet establish stable predictions
or repair the failed biological uncertainty coverage.
The subsequent [adaptive support fitter](Tools/Omics/CountObservation/Joint/Adaptive/README.md)
now bounds the continuous likelihood gap for **all 19 available models** from
two initial grids. All 19 also pass the 1% query-moment sensitivity criterion;
1.82 million independent numerical comparisons pass. This resolves the tested
support-convergence problem on the tested panel. The [full-gene extension](Tools/Omics/CountObservation/Joint/Adaptive/Full/README.md)
now processes all 15,706 Kang genes: 8,405 converge, nine eligible genes retain
solver limits, and 7,292 remain unavailable. Over 20.27 million independent
numerical comparisons pass; 44 Swift tests include a repaired real-count
convergence failure. HIRISA now processes all 18,082 genes: 13,082 converge,
185 retain the leaf-budget limit, and 4,815 are unavailable; 24.99 million
independent comparisons pass. Training-only dispersions are also verified for
all 13 donor omissions. Separate joint fits are now running; the first 64-gene
Kang and HIRISA pilot shards pass independent checks. These pilots do not
qualify the complete folds. The first two complete Kang donor omissions now
produce [control-only joint predictions](Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/README.md):
RMSE improves 20.7–30.8% versus no change and 21.2–25.6% versus training-mean
response on their eligible genes. These are reused development donors; the
parameter uncertainty and independent biological validation remain open; the
completed all-fold result below supersedes this initial two-fold snapshot.

The [complete 13-fold development evaluation](Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/evidence/2026-09-12-complete-reports/full-thirteen-fold-summary.json)
now includes all native fits, control-only predictions and independent numerical
and scoring checks. Kang improves RMSE 28.07% versus no change and 23.80% versus
training-mean response; all eight donors improve against both. HIRISA improves
50.84% versus no change but has **75.84% higher RMSE than training mean**, with
all five donors worse. The frozen pooled development criterion passes (29.05%
and 23.04% gains), while the HIRISA criterion fails. Of 216,058 source gene-folds,
131,342 are scored; 84,171 lack dispersion, 543 reach leaf limits and two do not
converge. These reused development donors do not establish independent biological
validation or calibrated uncertainty. The linked snapshot retains reports and
shard hashes; remaining raw evidence is retained locally pending archival capacity.

The [complete Replogle experiment](Tools/Omics/PerturbationPrediction/Replogle2020/RESULTS.md)
retains all 32,829 confident cells, 30 targets and 33,694 RNA features across five
technical gemgroups. All 150 native fits and predictions pass replay; all 750
output vectors agree with independent calculations. Inputs and predictions were
frozen before scoring, and repeated scoring is identical. [Gene identity evidence](Tools/Omics/PerturbationPrediction/Replogle2020/IDENTITIES.md)
uses a reproduced guide/gene table and reference sequence checks; full guide
sequences and genome-wide specificity remain unverified. The five groups share
K562/UPR biology and are not five independent biological replications.

The [HIRISA preparation-transfer test](Tools/Omics/Benchmarks/HIRISA/CONTEXT_TRANSFER.md)
now completes all 120 donor-excluded folds across 705,365 cells. Cross-preparation
ridge beats no-change in all twelve contrast means, but beats the simpler
cross-preparation mean in only **3/12**. All native replays and independent
numerical checks pass; this tests conditional RNA response between enriched
preparations and PBMCs in one study.

[Native Visium directory import](Tools/Omics/Multimodal/Visium/README.md) now
preserves the complete 4,039-spot lymph-node release directly from its original
counts and position table. All counts, identities and coordinates match Scanpy;
this closes a spatial input path without establishing tissue prediction.

[Native balanced reference mapping](Tools/Omics/ReferenceMapping/Logistic/README.md)
now reproduces the external classifier's labels for all 8,569 Baron query cells.
Macro-F1 improves over kNN in all four donors, while overall accuracy falls in
two and rare-class misses remain. These are candidate source labels with
uncalibrated probabilities; broader biological prediction claims remain unchanged.

The [complete Kang–Ding reference-transfer test](Tools/Omics/ReferenceMapping/CrossStudy/README.md)
now maps all 68,704 query cells through an explicit 14,976-gene panel while
preserving complete RNA-library denominators. Independent numerical checks pass,
but both coarse-family transfer targets fail (macro-F1 0.613/0.666), with zero
megakaryocyte recall in both directions. Cross-study annotation remains unqualified.
The [post-result source-label diagnosis](Tools/Omics/ReferenceMapping/CrossStudy/Diagnosis/README.md)
now verifies every original Kang label and all 68,704 cells' diagnostic scores.
The studies' same-named rare classes have sharply different RNA profiles; their
biological equivalence remains unqualified. Both original transfer failures stand.

The [single-cell program](Documentation/SingleCellInteroperability.md) now includes
native H5AD interchange, negative-binomial DE, sparse PCA/neighbors/clustering,
integration, marker scoring and multimodal count interchange. The complete
[1.61-million-cell HIRISA cohort](Tools/Omics/Benchmarks/HIRISA/README.md) has
published ingestion, DE, donor-response prediction, PCA, graph and seed-7
integration results. The original [program-preservation check](Tools/Omics/Benchmarks/HIRISA/INTEGRATION_PROGRAMS.md)
found decoder-dependent failures. A separately declared [within-library fit](Tools/Omics/Benchmarks/HIRISA/PROGRAM_CALIBRATION.md)
meets the same loss margins with 28/32 sensitive controls; four remain insufficient.
This development result preserves the original failures and does not establish
complete biological preservation. Full-cohort
[clustering publication, replay and independent checks](Tools/Omics/Benchmarks/HIRISA/FULL_CLUSTERING_RESULTS.md)
now pass. The original 9.1-hour clustering baseline has also finished with an
identical complete result; its separate reconstruction failed for disk space.
[Native program scoring](Tools/Omics/Benchmarks/HIRISA/NATIVE_PROGRAM_RESULTS.md)
also reproduces the complete cohort’s reference scores exactly and passes replay.
The new [annotation-retention diagnostic](Tools/Omics/Benchmarks/HIRISA/ANNOTATION_RETENTION.md)
finds native losses beyond its margins in **37/146** supported, sensitive
comparisons, including **10/29 rare-label comparisons**. All three Harmony
references also fail the complete gate. Broader biological preservation remains
unqualified. A single [protected-stratum regression experiment](Tools/Omics/Benchmarks/HIRISA/PROTECTED_REGRESSION.md)
retains 37/146 failures: four resolve and four new failures appear. This candidate
has not been promoted into the native solver.
A separate [exact-tree MNN trial](Tools/Omics/Reduction/MNN_TREE.md) preserved all
outputs on three original cohorts but ran slower, so its prototype remains archived.
A [local-kernel MNN candidate](Tools/Omics/Reduction/LocalMNN/README.md) is faster
on the two larger cohorts but loses cell-type and program signal, so it also
remains experimental.
The [tiled all-anchor Gaussian option](Tools/Omics/Reduction/GaussianKernel/README.md)
runs about **2.4× faster** than the preserved production scalar owner on all
three cohorts, with identical measured biological metrics and neighbor/prediction
arrays. The scalar default is also faster and retains exact original outputs.
These are native-owner timings; independent-validation limits remain.
A fixed [approximate-matching/full-Gaussian follow-up](Tools/Omics/Reduction/FullGaussianMNN/README.md)
restores the measured Ding gates but still fails Kang megakaryocyte recall.
The approximate matcher remains experimental; production retains exact matching.

## Input preservation and execution evidence

Native [H5AD axis projection](Tools/Omics/H5AD/Projection/README.md) now reads the
original legacy Kang file directly: all 24,673 cells × 15,706 genes, annotations,
category definitions and PCA/UMAP values agree with AnnData. Full and repeated
selections reconstruct exactly; earlier unique and repeated projections remain
byte-identical. This closes a specific input-preservation gap without changing
the biological prediction evidence above. Other legacy formats and analytical
identity restrictions remain explicit.

The [complete original-data analytical route](Tools/Omics/H5AD/Projection/LEGACY_COUNT_ROUTE.md)
also passes: native annotation, all 14,184,532 count records and all 124
pseudobulk groups agree with the original Kang information. RNA assay totals
match every cell, and donor/treatment identities remain explicit. This verifies
input preservation and arithmetic, without adding a new prediction claim.

The [full HIRISA annotation route](Tools/Omics/H5AD/Annotation/Atlas/README.md)
now preserves the 6.1 GB atlas while adding source-row provenance to all
1,612,594 cells. All 46 original datasets and 7.74 billion stored elements match;
the native run used 180.9 MiB peak resident memory on an Apple M4. APFS publication
and shallow group detachment avoid duplicating untouched data. This qualifies
annotation/storage handling, not new biological labels or predictions.

A [real-Hagai Metal kNN distance-block experiment](Tools/Omics/Reduction/MetalDistanceBlocks/README.md)
now reduces research-harness median elapsed time by 35.1% versus scalar FP32 CPU
on all 13,863 cells. Neighbor membership matches the full FP64 reference, with
one order difference and small distance errors. This is not yet a product
backend, full-pipeline speedup or biological-preservation qualification.

The first [Metal count-normalization check](Tools/Omics/CountStore/Metal/README.md)
now covers every one of the original Kang dataset's 14,184,532 records on physical
M4/M4 Pro GPUs. The explicit FP32 option passes its declared numerical tolerance
and native replay; the FP64 CPU default remains byte-exact. End-to-end medians
are 1.043 s CPU and 1.050 s Metal, so this experiment establishes **no speedup**.
A subsequent [shared-writer qualification](Tools/Omics/CountStore/Metal/Profile/README.md)
preserves every output byte and lowers final medians to 0.918 s CPU / 0.912 s
Metal. The first CPU call is slower and its three-run mean is slightly worse;
these are bounded measurements, not a general speedup claim. The real single-cell
CLI and bounded memory-safety checks pass; downstream biological claims are unchanged.

## One scientific question, several scales

NumiVivo is being developed to connect a molecular hypothesis to a dynamic biological model:

```text
Structures and force fields
        ↓
Molecular dynamics and sampled configurations
        ↓
Electronic structure, QM/MM and orbital embedding
        ↓
Reaction energetics and explicitly qualified kinetic evidence
        ↓
Reaction networks, target engagement and physiological dynamics
        ↓
Versioned experiments, observations and reproducible results
```

This is the integration direction, not a claim that every arrow is already an automated, validated calculation. In particular, an electronic-energy difference is not automatically an activation free energy, reaction rate, binding affinity or biological outcome.

## Explore the capabilities

| Area | What the source provides | Explore |
|---|---|---|
| **Molecular foundations** | Canonical atoms, bonds, residues, conformers and periodic cells; selections and atom mapping; PDB, mmCIF, SDF/MOL V2000, MOL2 and a strict SMILES subset; topology preparation and unit-explicit force-field records. | [Structure and force-field foundation](Documentation/Design/MOLECULAR_FOUNDATION_WAVE_A.md) |
| **Apple-native molecular dynamics** | Metal force kernels, neighbor construction, minimization, NVE and Langevin NVT, molecular-center Monte Carlo NPT, PME electrostatics, distance constraints and linear virtual sites. Preparation protocols retain stage identity, restart state and bounded trajectory output. | [MD protocols](Documentation/Design/MD_PROTOCOL_WORKFLOW.md) |
| **Native electronic structure** | FP64 Cartesian Gaussian integrals, restricted/unrestricted Hartree–Fock, restricted LDA, embedded Hamiltonians, MP2 and small-system configuration interaction. The chemistry path runs without Python callbacks or a CUDA runtime. | [Native chemistry example](Examples/native-chemistry/README.md) |
| **Embedding and reaction research** | QM/MM electrostatic and boundary-link machinery, C-PCM reaction-field work, correlated orbital information, orbital-subspace alignment and a bounded path-consistent QIO optimizer. These are experimental methods, not a reproduced protein-reaction result. | [Embedding source](Sources/NumiVivoKit/Embedding) · [QM environment source](Sources/NumiVivoKit/QMEnv) |
| **Programmable reaction dynamics** | Typed molecular programs and compiled ProgramPacks; deterministic kinetics; discrete stochastic simulation; exact-SSA/tau-leap/RK2 execution across dependency-separated components; temporal rules, monitors and concentration transport within declared backend limits. | [Hybrid runtime](Examples/hybrid-reaction-runtime/README.md) · [ProgramPack backend](Documentation/Design/PROGRAM_PACK_METAL_BACKEND.md) |
| **Target engagement and physiology** | Exposure-driven reversible binding, covalent conversion, competition and turnover; a native FP64 reference and an existing-runtime Metal path; physiological exchange and molecular–physiology coupling contracts. | [Target-engagement example](Examples/target-engagement/README.md) |
| **Single-cell analysis and prediction** | Native H5AD/H5MU count interchange, paired negative-binomial DE, sparse reduction and integration, marker programs, donor-response and target/composition baselines. | [Single-cell workflows](Documentation/SingleCellInteroperability.md) · [Prediction evidence and limits](Documentation/BiologicalPrediction.md) |
| **Reproducible experiments** | Content-addressed artifacts and tasks, evidence references, configuration-bound checkpoints, staged protocols, compact trajectory chunks, observation records and explicit failure results. | [Artifacts and provenance](Documentation/Design/ARTIFACTS_AND_PROVENANCE.md) · [Trajectory storage](Documentation/Design/MD_TRAJECTORY_ARCHIVE.md) |

The repository also contains experimental quantum-algorithm, reaction-path, reaction-network, population, calibration and surrogate components. Their presence is not a claim of complete Qiskit, ORCA, GROMACS or physiological-modeling equivalence. See the [source-backed capability map](Documentation/CAPABILITIES.md) before choosing a backend.

## Prepared systems and chemistry-to-rate workflows

The [prepared molecular workflow guide](Documentation/PreparedMolecularWorkflows.md) connects explicit protonation/stereo preparation and native parameter assignment, checkpointed adaptive replica sampling, density-fitted active spaces, mapped reaction connectivity and context-bound conditional rates through the existing runtimes and artifact DAG. It includes the new `molecule-*` commands and `h3-connected-rate` workflow fixture. Accepted replica states can be [exported and freshly verified](Documentation/Design/MOLECULAR_SAMPLING_EXPORT.md) as typed workflow inputs, retaining exact checkpoints, atom mapping and partial sampling status. The [native/CLI audit](Documentation/Audit/2026-09-08_SAMPLING_EXPORT.md) records the tested scope.

These integrations do not yet provide arbitrary-triclinic Metal MD, automatic pKa/CIP chemistry, a universal native force field or complete protein reaction free energies. The guide distinguishes implemented source paths from the remaining physical and execution qualification.

The [prepared molecule to observable example](Examples/prepared-reaction/README.md) now supplies a complete finite route: prepare a synthetic harmonic H₂ system, sample an accepted state, use its geometry for fresh H₃ exchange qualification, and calculate a tagged atom's conditional reaction probability. Its [native and public-CLI audit](Documentation/Audit/2026-09-08_PREPARED_REACTION_CAMPAIGN.md) retains the model transfer, partial sampling status, reaction evidence and explicit bath assumptions.

The [independent MD benchmark panel](Tools/Benchmarks/README.md) now checks realistic
water, protein, ligand, DNA and membrane systems. Its [audit](Documentation/Audit/2026-09-08_FRONTIER_MD_REFERENCE_PANEL.md)
records corrected torsion forces, opt-in compensated positions, RATTLE integration,
and short energy-conservation/refinement results. Failed inputs and the remaining
equilibrium, performance and broader-model qualifications remain explicit.

## Built around Apple silicon

**Native execution rather than a Python simulation loop.** Swift manages scientific objects, experiments and concurrent operations. C++23 implements compilation, validation and portable reference components. Metal executes molecular-dynamics and kinetic kernels. Precision-sensitive chemistry retains FP64 native CPU calculations rather than forcing every calculation onto the GPU.

**State stays where the calculation needs it.** The GPU runtimes use persistent buffers, compiled tables and explicit capacities. Private simulation state is separated from host-visible commands, diagnostics and requested observations. Full coordinate readback is an explicit sampling or checkpoint operation.

**A failed candidate must not become accepted state.** Runtime transactions distinguish proposed and accepted evolution. Checkpoints bind the represented state to its model and numerical configuration. The MD numerical profile is versioned separately so a restart cannot silently substitute a different algorithm.

**Reproducibility has a defined scope.** Counter-based random namespaces retain seeds and accepted-step identity. They do not guarantee bitwise-identical trajectories across devices or compiler versions; PME accumulation and floating-point execution have additional reproducibility limits. Integrity hashes establish which bytes were used, not whether the scientific model is correct.

## Get started

The primary development target is an **Apple-silicon Mac with macOS 15 or newer**, Swift 6 and an Apple SDK providing Metal and a C++23-capable toolchain. The package also declares iOS 18 support; it is not a released iOS application. The native package has no Python or CUDA runtime dependency.

Run from a terminal on the Apple machine:

```sh
git clone https://github.com/Numi2/numiVivo.git
cd numiVivo
swift build -c debug
.build/debug/numivivo --help
```

These are build and execution instructions, not a recorded successful build of the current revision. Keep the commit SHA with compiler diagnostics and numerical results. Portable checks do not qualify the complete Apple package.

### First experiment: exposure and target occupancy

Use a supplied synthetic fixture—no downloaded dataset or prepared protein is required:

```sh
RUN_DIR="$(mktemp -d "${TMPDIR:-/tmp}/numivivo.XXXXXX")"

.build/debug/numivivo engagement-validate \
  Examples/target-engagement/synthetic-pulse.json

.build/debug/numivivo engagement-run \
  Examples/target-engagement/synthetic-pulse.json \
  --backend reference \
  --output "$RUN_DIR/engagement.json"
```

Inspect target-state fractions over the supplied exposure schedule. The inputs and observations are **synthetic mathematical fixtures**, not measured pharmacology or treatment recommendations. The [example guide](Examples/target-engagement/README.md) adds the Metal backend, compiled kinetic source, conditional free-energy conversion and held-out observation evaluation.

## Choose an experiment

| Start with | What it exercises | Input |
|---|---|---|
| [Native chemistry](Examples/native-chemistry/README.md) | Gaussian integrals → Hartree–Fock → embedded Hamiltonian → MP2/FCI → orbital information. | Included H₂/STO-3G example. |
| [Prepared molecule to observable](Examples/prepared-reaction/README.md) | Preparation → Metal sampling → verified geometry seeds → fresh reaction qualification → conditional probability. | Published synthetic H₂ library, H₃ exchange seeds and maintained H₂ bath. |
| [Hybrid reactions](Examples/hybrid-reaction-runtime/README.md) | Exact SSA, tau-leaping and RK2 on separate reaction components; checkpoint and resume. | Included synthetic model and counts. |
| [Target engagement](Examples/target-engagement/README.md) | Exposure, binding, covalent conversion, competition and turnover. | Included synthetic experiments. |
| [Prepared-system MD](Examples/md-preparation/README.md) | AMBER import → minimization → NVT → NPT → production samples. | Your prepared topology and restart. |
| [Digital tissue homeostasis](Examples/digital-tissue-homeostasis/README.md) | Molecular-control, host-context and coupling concepts. | Included research fixtures; check the declared runtime support. |

### A prepared molecular system, one MD protocol

For an appropriately prepared **orthogonal periodic system**, the existing AMBER bridge provides the starting point. Once `system.json` and `initial.json` have been imported:

```sh
.build/debug/numivivo md-protocol-template system.json > protocol.json

# Review the durations, temperature, pressure and minimization gate first.
.build/debug/numivivo md-protocol-validate protocol.json \
  --system system.json --state initial.json

.build/debug/numivivo md-protocol-run protocol.json \
  --system system.json --state initial.json \
  --store ./md-artifacts > receipt.json
```

The receipt identifies durable restart state. Coordinate samples are stored as bounded binary chunks, not an ever-growing in-memory JSON trajectory. Later stages preserve state unless their velocity initialization is explicitly changed. A failed required minimization or rejected dynamics candidate blocks the protocol rather than silently changing the experiment.

The template is illustrative: its duration does not establish equilibration. [Import, sampling and resume details →](Documentation/Design/MD_PROTOCOL_WORKFLOW.md)

## Molecular biology inside NumiLab

The wider goal is a coupled experimental environment in which molecular reactions respond to cell state, tissue organization, physical transport and nervous-system activity.

| System | Intended responsibility |
|---|---|
| **NumiVivo** | Molecular representations, chemistry, reaction dynamics and molecular experiments. |
| **NumiTissue** | Cell populations, development and tissue organization. |
| **NumanX** | Geometry, mechanics, fluids and physical coupling. |
| **NumiBrain** | Neural, autonomic and embodied control. |

NumiVivo contains coupling contracts and participant infrastructure. Fully qualified cross-repository biological simulations remain an integration objective, not something established by this table.

## Development direction

The immediate development priority is **the experimentally evaluated single-cell
prediction workflow**: H5AD → real-data benchmarks → native negative-binomial DE
→ PCA/neighbors/clustering → integration with biological preservation →
perturbation prediction. Finish full-cohort qualification, independent-study
prediction and context-transfer evidence before claiming general biological
prediction. See the [current assessment and next gates](Documentation/BiologicalPrediction.md#next-evidence-needed).

The broader molecular program remains an integrated, numerically qualified
research workflow.

The priorities are to consolidate model semantics and execution ownership; qualify the Apple MD and kinetic backends; expand electronic-structure and embedding methods against independent references; and connect reaction energetics to observations and kinetics without losing thermodynamic meaning. Membrane interfaces, general triclinic dynamics, delayed/refractory state, live fidelity migration and larger correlated calculations require additional work.

The chemistry program is informed by [CovAngelo, arXiv:2604.10487](https://arxiv.org/abs/2604.10487). Its methods motivate native QM/MM, compact embedding and consistent reaction-path treatment. This repository does not claim to reproduce that paper's results or contain its authors' unpublished implementation.

## Documentation and contribution

[Documentation index](Documentation/README.md) · [Capability map and limits](Documentation/CAPABILITIES.md) · [Execution audit](AUDIT.md) · [Contribution guide](CONTRIBUTING.md) · [Security policy](SECURITY.md)

```text
Sources/NumiVivoCore/      C++ compiler, validation, pack format and reference logic
Sources/NumiVivoKit/       Swift scientific modules, runtimes, artifacts and workflows
Sources/NumiVivoShaders/   Metal kernels and explicit shader-module loading
Sources/NumiVivoCLI/       numivivo command-line interface
Examples/                 Model fixtures and executable examples
Schemas/                  Versioned data contracts
Documentation/            Methods, architecture, audits and implementation limits
Tools/                    Development tools and portable qualification harnesses
```

Report reproducible bugs through the repository's Issues tab, including the commit, command, platform and relevant diagnostics. Keep confidential research data out of public reports. Contributions should include the affected contract, an example and an explicit validation boundary; see [CONTRIBUTING.md](CONTRIBUTING.md).

## Scientific boundary and license

NumiVivo models biological systems; it does not validate a therapy, authorize an experiment or establish safety in a living organism. An executable model is not proof of a realizable biological construct. Numerical checks, calibration evidence and biological validation are distinct.

Licensed under [Apache License 2.0](LICENSE). See [NOTICE](NOTICE). External datasets, force-field parameters, model weights and third-party materials retain their own terms. For research use, cite the exact repository revision and the methods and source data used by the calculation.

The completed donor-exclusion evaluation now also has a [post-hoc support diagnostic](Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/evidence/2026-09-12-support-diagnostic/README.md): HIRISA’s loss to the training-mean baseline persists within the training rate range in every donor. General biological outcome prediction remains unestablished.

A [response-transport development experiment](Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/ResponseTransport/README.md) is now complete: it reduces HIRISA error but still fails against the mean-response baseline and worsens Kang against the original joint model. It is not promoted; independent biological validation remains open.

The [Parse prediction endpoint audit](Tools/Omics/PerturbationPrediction/ParseIFNB/BCellCounts/NEXT_STEPS.md#prediction-endpoint-audit-2026-09-12) identifies a remaining input gap: validated pseudobulk count sums cannot substitute for the frozen model’s mean per-cell log-normalized expression. All twelve donor reports were checked; no Parse prediction has been scored.

The [streaming log-CPM accumulator](Tools/Omics/PerturbationPrediction/StreamedLogCPM/README.md) now passes complete real-matrix numerical verification (maximum mean error 1.95e-14). The complete scoped library build and linked stream harness also pass; the [twelve-donor Parse replay](Tools/Omics/PerturbationPrediction/StreamedLogCPM/ParseReplay/README.md) and full-cohort terminal review pass for all 72,446 B cells and 40,352 features (maximum error 7.99e-15). No prediction claim follows from these numerical checks.

The [native feature-major log-CPM stream](Tools/Omics/PerturbationPrediction/StreamedLogCPM/FeatureMajor/README.md) now handles all 361.6 million Norman records without reordering. All 7.99 million condition means pass independent verification; large result arrays remain resident.

The [feature-at-a-time log-CPM bundle](Tools/Omics/PerturbationPrediction/StreamedLogCPM/FeatureStream/README.md) now processes the complete Norman matrix at 44.1 MB peak CLI RSS without resident condition-by-feature results. All means and axes are exact; the lower-memory path is slower in the recorded run.
