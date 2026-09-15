# NumiVivo completion roadmap

The product objective is a leading Apple-native tool for reproducible molecular-to-biological research. The scope is the connected workflow described in the README, capability map and engineering principles. Passing one build, one molecular fixture or one numerical suite does not complete that objective.

This roadmap orders the work without replacing the intended architecture with a smaller product. Unsupported methods must remain explicit, and scientific claims must retain their actual evidence class. An unavailable author dataset can block a particular reproduction; it does not block improving the general platform.

## Acceptance requirements

| Requirement | Evidence needed for completion | Current work or remaining gap |
| --- | --- | --- |
| Reliable execution and recovery | Clean full Apple build; verified output/cache identities; cancellation and concurrent caller tests; injected storage-failure recovery; uninterrupted versus resumed results across complete protocols. | September 2026 work closes shared-task cancellation and MD resume/publication defects. Native verification is recorded separately; these fixes do not qualify every runtime. |
| Correct molecular and stochastic dynamics | Independent force/energy checks, conservation and constrained-ensemble checks, PME convergence, exact stochastic distributions, hybrid authority conservation, transactional rejection, and restart invariance over the supported numerical profiles. | Seven dedicated hybrid native conformance tests passed on M4 Pro; their exact scope is recorded in the reliability audit. The [v3 numerical audit](Audit/2026-09-08_MD_NUMERICAL_QUALIFICATION.md) adds independent charged-PME, finite-time Langevin, fixed-geometry constrained thermalization and ideal molecular NPT transition gates. The [v6 interacting audit](Audit/2026-09-08_INTERACTING_ENSEMBLES.md) adds seven-system 5 ps conservation/refinement and a bounded rigid-water NVT pass at 0.5 fs over 210 ps. Broader ensembles, interacting NPT and force-provider combinations remain open. |
| Qualified chemistry and sampled uncertainty | Independent electronic references, mapped path/connectivity evidence, declared free-energy and nuclear approximations, convergence studies, surrogate held-out errors and exact acceptance accounting; uncertainties remain correlated and missing sources remain visible. | The reaction gate is repaired at `9f9da42`. Revision `59e1868` adds bead-convergence, correlated-sampling and surrogate-coverage assessments; their declared scope is in the nuclear/surrogate contract. No protein-scale kinetic or general chemical accuracy claim follows from small-system gates. |
| Connected research workflows | Complete prepared-system → sampling → electronic/path calculation → conditional kinetic evidence → observable workflow, with units, mapping, context, parameter provenance and assumptions retained at every edge; both successful and missing-input cases. | The [prepared-reaction milestone](Audit/2026-09-08_PREPARED_REACTION_CAMPAIGN.md) completes a published finite hydrogen route through both native SDK and public CLI, including missing-input and altered-evidence failures. Its generic handoff interfaces retain exact accepted geometry, fresh reaction qualification and conditional reservoir assumptions. Larger, realistic and ensemble-averaged research claims still require their own validation. |
| Complete declared runtime semantics | Source, native tests and counterexamples for the accepted language/backend contract; persistent state includes everything that can affect continuation; coupled participants commit the same logical boundary. | ProgramPack delays/refractory behavior, broader spatial interfaces, live authority migration and stronger multi-participant execution remain separate from the well-mixed hybrid backend. Follow the linked design contracts rather than interpreting fidelity names as implementation evidence. |
| Demonstrated useful scale and performance | Published systems, source/binary/device identities, correctness tolerances, preparation and warmup policy, elapsed time, throughput, memory, and accepted-work accounting; compare like-for-like numerical models against named independent references where available. | No broad speedup or maximum biological scale is established. Measure complete workflows and bottlenecks before promoting an acceleration claim. |
| A usable, maintainable research product | Documented complete examples, discoverable CLI/API operations, actionable bounded errors, readable progress and result inspection, automated native regression gates, install/release checks, and retained evidence for failures as well as successes. | Workflow activity now has a concrete ownership model. CLI examples, operation coverage, long-run inspection and release qualification must stay aligned with the actual package. |

## Immediate single-cell execution order

The active development sequence is H5AD/AnnData → complete real-data benchmarks
→ native negative-binomial DE → PCA/neighbors/clustering → integration with
biological preservation → perturbation prediction. The
[twelve-part single-cell roadmap](SingleCellInteroperability.md) retains all
required capabilities, including annotation, multimodal assays, out-of-core
execution, later Metal acceleration and cross-scale biology.

The benchmark stage now has a checked-in [scope and provenance registry](../Tools/Omics/Benchmarks/benchmark-manifest.json). Its offline checker validates the seven declared public/technical entries and retains the current conclusion that biological-outcome prediction is not established.

The [million-cell AnnData annotation and graph-label export](../Tools/Omics/H5AD/Annotation/Atlas/Clustering/README.md)
now complete with exact cell-identity/metadata checks and backed AnnData reopening.
Numeric legacy categories, scalar embeddings and numeric-index analytical
handoff are also implemented with their declared conformance checks; remaining
format and consumer bounds are listed in the interoperability roadmap. These
completed repairs should not be scheduled again as missing capabilities. The
external gzip-wrapped H5AD path now has native import, streamed scan, axis
projection and annotation qualifications. The [annotation receipt](../Tools/Omics/H5AD/Import/Gzip/evidence/2026-09-14-annotation/README.md)
binds the exact compressed source and publishes a new plain H5AD after bounded
private decoding; this remains interoperability/storage evidence, not biological
annotation or outcome validation.
The [opt-in Metal sparse feature-statistics profile](../Tools/Omics/Reduction/MetalDistanceBlocks/NativeOwner/evidence/2026-09-14-feature-statistics/README.md)
now feeds HVG/PCA moments from sparse nonzero values on a physical Apple GPU,
with CPU FP64 still the default. It is a bounded numerical qualification and
does not close general Metal acceleration or biological-outcome validation.
The [opt-in Metal sparse PCA-operator profile](../Tools/Omics/Reduction/MetalDistanceBlocks/NativeOwner/evidence/2026-09-14-pca-operators/README.md)
now runs resident Krylov projection and transpose over sparse CSR/CSC streams on
a physical Apple GPU. Its three focused tests pass with an explicit FP32
residual tolerance; CPU FP64 remains the default and broad acceleration remains
unqualified.
The [opt-in Metal batched NB objective profile](../Tools/Omics/Reduction/MetalDistanceBlocks/NativeOwner/evidence/2026-09-14-nb-objective/README.md)
now evaluates the mean-dependent NB2 line-search objective in bounded FP32
Metal batches. Count-only gamma terms, coefficient updates, diagnostics and
the published log likelihood remain on the exact FP64 CPU owner. Its four
focused native tests pass, but this is only an objective component: full GPU
model fitting, timing, FDR/interval calibration, million-cell scaling and
biological-outcome validation remain open.

The current [PBMC3K reduction receipt](../Tools/Omics/Reduction/MetalDistanceBlocks/EndToEnd/evidence/2026-09-15-pbmc3k-current/README.md)
now closes the active PCA → neighbors → clustering owner check on a pinned
2,700-cell public source. Native CPU/Metal and independent Scanpy agree on the
selected feature IDs and all 77,916 final graph edges; binary graph verification
and three seeded Louvain replays agree across native backends. A real feature-ID
namespace failure is retained with its explicit `gene_ids` repair. This is
numerical interoperability and descriptive clustering for one unreported
library, not biological-outcome validation.
The current [Kang donor-integration receipt](../Tools/Omics/Reduction/evidence/2026-09-15-kang-donor-integration/README.md)
also completes the fitted-PCA donor correction and independent residual checks
on eight donors. It remains transductive and does not establish held-out donor
prediction or preservation of biological signals.

The [AlphaGenome Atlas boundary](AlphaGenomeAtlas.md) remains hypothesis-only:
its variant/regulatory molecular-effect scores do not provide measured RNA,
protein, phenotype or tissue outcomes. The native cross-scale contract therefore
still has no complete qualified path, even though the bounded VCF/Atlas
provenance adapter and offline checks are implemented.

The [prediction assessment](BiologicalPrediction.md) records the current
scientific answer: bounded expression-response prediction is demonstrated,
while reliable general biological-outcome prediction remains open. Complete
HIRISA ingestion, DE, 79 donor folds, PCA/graph and seed-7 integration have
published evidence. [Full-cohort clustering](../Tools/Omics/Benchmarks/HIRISA/FULL_CLUSTERING_RESULTS.md)
now completes native publication, replay and independent partition checks on all
1,612,594 cells. This qualifies graph partitions, not calibrated cell types.
The [preparation-transfer experiment](../Tools/Omics/Benchmarks/HIRISA/CONTEXT_TRANSFER.md)
also completes all 120 frozen folds: cross-preparation ridge beats both required
baselines in only 3/12 contrast means and is worse than matched within-preparation
ridge in all twelve. These completed experiments are no longer pending work.

The thirteen [donor-excluded joint count-response folds](../Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/README.md)
are now fully fitted, predicted, independently verified and scored under the
frozen endpoint plan. Kang improves RMSE 23.80% over training mean across eight
donors; HIRISA is 75.84% worse and loses in all five donors. The pooled criterion
passes, but the HIRISA criterion fails. These reused development donors do not
provide fresh external validation or repair the failed historical interval coverage.

Training-only response diagnostics, donor-selected shrinkage, native two-study
transfer and the range/error decomposition are complete. The
[three-study held-out follow-up](../Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/ResponseShrinkage/StudyHeldOut/README.md)
now selects penalties by held-out training studies with equal study weights.
All 75 donors and 885,000 predictions pass independent numerical verification.
It still fails HIRISA in every donor; no study reaches a 5% gain over both
baselines. These are completed development tests on previously inspected
cohorts, including explicit reuse of GSE181897. They are not pending jobs or
fresh external qualification.

The [multigene context-kernel follow-up](../Tools/Omics/PerturbationPrediction/StudyContextKernel/README.md)
now uses the full control profile to predict shared response weights. Native
execution and independent verification complete for all 75 donors and 885,000
predictions. HIRISA improves 21.11% over the training mean, but Kang and
GSE181897 lose to that baseline. The original gene-wise failures remain retained;
neither model establishes across-study reliability. Further work must represent
relevant study context and uncertainty, with a new external protocol frozen
before outcome inspection. Do not select a winner per inspected study or change
the required baselines. Full raw joint-shard archival remains incomplete; the
context-kernel and shrinkage experiments retain complete declared archives.

The [Kang-HIRISA NB2 follow-up](../Tools/Omics/PerturbationPrediction/CrossStudyIFNB/NB2_CROSS_STUDY.md)
now runs the opt-in paired-donor count model on all 26 cross and within-study
folds. Native execution and independent source-count scoring pass all 130
estimate vectors. HIRISA → Kang gains 2.98% over no change but loses to the
cross-study mean; Kang → HIRISA is 0.56% worse than no change. The result
strengthens the implementation evidence without closing context transfer or
general biological-outcome prediction.

Integration still requires resolving sensitive marker/program and rare-cell
preservation losses. The [full-HIRISA rigid projection](../Tools/Omics/Benchmarks/HIRISA/RigidProjection/README.md) now tests all original cells with unchanged diagnostics: 31/146 annotation failures remain, including 9/29 rare failures. Within-donor isometry does not repair cross-donor preservation; the method is not promoted. The [source design audit](../Tools/Omics/Reduction/IntegrationDesign/README.md)
now verifies full-atlas donor/batch crossing, subset rank deficiencies, and exact
monocyte-preparation/batch overlap. The [protected sample-group admission guard](../Tools/Omics/Reduction/IntegrationProtection/README.md) now rejects disconnected
joint condition/group designs in ridge and MNN. An opt-in additive donor+batch
ridge path also checks pairwise factor co-occurrence and retains factor-major
axes, but it is bounded to two covariates and remains transductive numerical
development evidence. It does not constrain correction or establish biological
preservation, interactions or unseen-donor mapping; retain the existing
scientific gates. The primary Adamson paper and Table S1 now confirm the two UPR controls and their sequences. A separate [provisional role-sensitivity run](AdamsonProvisionalRoleSensitivity.md) fits and scores 80 supported held targets under an explicitly pooled `pBA580`/`pBA582` control assumption, but deposited-label mapping and the 94-versus-93 roster discrepancy remain open before authoritative fitting of its frozen transfer test, and Parse's intervention and feature identities before admitting its proposed prediction test. Retain all negative results and simple baselines. Measured phenotype links need their own acceptance experiments.

The [condition-stratified donor candidate](../Tools/Omics/Benchmarks/HIRISA/CONDITION_STRATIFIED.md)
was executed on the complete HIRISA bundle and passes native replay plus an
independent per-cell SVD oracle. Its absolute retention result still fails
40/146 supported comparisons, including 10/29 rare comparisons, versus 37/146
for the existing native correction. The candidate remains a rejected development
experiment; local composition/alignment and independent biological endpoints are
still required.

The predeclared [ridge-4 global candidate](../Tools/Omics/Benchmarks/HIRISA/RIDGE4.md)
then changed only the fixed donor-effect penalty from 1 to 4. Native replay and
the independent SVD oracle passed, but absolute retention remained 37/146
failures, including 10/29 rare comparisons, with a slightly worse worst-fold
loss. The preservation gate remains open and no correction method is promoted.

The following molecular milestones remain part of the wider product roadmap;
they do not replace this immediate single-cell priority.

## Molecular qualification sequence

The [frontier benchmark contract](Design/FRONTIER_BENCHMARKS.md) and
[reference campaign tools](../Tools/Benchmarks/README.md) make the next broad
simulation milestone executable: realistic-system force/energy comparisons,
explicit unsupported inputs and short native execution checks. These gates
precede interacting ensemble, triclinic, performance, sampling and broader
chemistry qualification; a completed scorecard does not close those milestones.

The [implemented frontier MD milestone](Audit/2026-09-08_FRONTIER_MD_REFERENCE_PANEL.md)
now includes the independent panel, torsion-force repair, compensated coordinates,
RATTLE correction and short smooth-model NVE conservation/refinement. The
[longer interacting follow-up](Audit/2026-09-08_INTERACTING_ENSEMBLES.md) passes
5 ps conservation/refinement on seven prepared systems and the prescribed global
kinetic/configurational checks for rigid water at 0.5 fs over 210 ps. The original
1 fs kinetic failures and failed water preparation remain retained. Broader NVT,
NPT, cell-shape, performance, sampling, chemistry and release gates stay open.

1. Completed the cancellation, MD recovery and hybrid conformance milestone with [native Mac mini evidence](Audit/2026-09-07_COMPLETION_RELIABILITY.md): a clean complete build, 224 baseline tests, 22 final-source regressions and 57 CLI checks. This milestone does not complete the acceptance table.
2. Completed the [streaming archive recovery milestone](Audit/2026-09-07_ARCHIVE_STREAMING.md): a real 100,001-chunk archive validates, resumes and extends through production code while preserving its prefix. Fresh CLI processes verify the archive with memory independent of its length for the measured fixture. The explicit materializing index retains its bound. Turn new integrity or recovery failures into owning-runtime repairs; do not relax gates to obtain a pass.
3. Extend the [v3 MD numerical milestone](Audit/2026-09-08_MD_NUMERICAL_QUALIFICATION.md) into the wider force/ensemble/continuation matrix. Nuclear/surrogate and resource-admission work remains evidence-scoped; respect physical GPU ownership.
4. Completed the representative [prepared, mapped full-route campaign](Audit/2026-09-08_PREPARED_REACTION_CAMPAIGN.md) using published inputs. Extend its scientific coverage while retaining source/destination models, mass/temperature differences, molecularity, standard states and transfer assumptions. The bimolecular encounter interface remains separate from the bound-complex unimolecular rate bridge.
5. Use measured correctness and resource bottlenecks to prioritize scale improvements and broader runtime semantics. Keep permanent interfaces and one numerical authority per represented state.
6. Audit the entire acceptance table against current source, native runs, artifacts and release behavior before claiming product completion.

## Source contracts

- [Capability map](CAPABILITIES.md) and [engineering principles](PRINCIPLES.md).
- [General workflows](Design/GENERAL_WORKFLOW_PLATFORM.md), [MD protocols](Design/MD_PROTOCOL_WORKFLOW.md), and [trajectory archives](Design/MD_TRAJECTORY_ARCHIVE.md).
- [Hybrid execution](Design/EXECUTABLE_HYBRID_RUNTIME.md), [ProgramPack backend](Design/PROGRAM_PACK_METAL_BACKEND.md), and [transactional runtime](Design/TRANSACTIONAL_RUNTIME.md).
- [Prepared molecular workflows](PreparedMolecularWorkflows.md), [QM/MM free energy](QMMMFreeEnergy.md), and [nuclear, uncertainty and surrogate methods](Design/NUCLEAR_UNCERTAINTY_SURROGATE.md).

Measured validation records must identify the tested source and scope. Historical `VALIDATION_STATUS.md` results remain evidence for their named revisions, not an automatic certificate for later code.
