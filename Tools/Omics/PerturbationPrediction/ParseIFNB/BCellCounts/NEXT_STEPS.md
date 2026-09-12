# From validated B-cell counts to a predictive test

The count workflow and complete archive review passed, published in c4b4f1dd.
This establishes faithful selected count processing; it does not establish
predictive accuracy. Preserve the frozen execution and its evidence.

## Final count evidence review

Before publishing a completed count result, require all of the following:

1. Twelve distinct frozen donors in each phase, with 72,446 cells and
   124,909,573 records per phase; no absent donor or failed attempt substituted.
2. Matching preparation, plan, runtime and source identities. Recheck the actual
   files against their recorded hashes, rather than trusting completion flags.
3. Revalidation of every native bundle against the independent cell totals and
   donor/condition feature aggregates. Require native exit success and the
   expected stream byte count for each donor.
4. Matching ingest/replay stream hashes and every fetched source range hash.
   Retain both phases, attempt logs and all source QC discrepancies.
5. An archive inventory with file sizes and SHA256 values, followed by a separate
   verification of the archive. Identify dependencies not embedded in it.

If any check is missing, publish the achieved partial state and its specific gap.
The terminal JSON flags alone are insufficient evidence. Preserve failed attempts
even when a later attempt succeeds.

## Admission before any prediction score

- The primary supplement identifies IFN-beta at 100 ng/mL, Biotechne / R&D
  Systems catalogue 8499-IF-010/CF; see [dose provenance](dose-provenance.json).
  Resolve cross-study protocol comparability separately. Do not infer an IU
  conversion or equivalent biological exposure from mass concentration.
- Preserve the versioned 11,600-feature exact-match contract and its renewed
  development results: HIRISA passes, while Kang and GSE181897 fail the transfer
  gate. This resolves source compatibility but does not promote the candidate.
  Do not silently zero-fill, alias or choose a panel after inspecting outcomes.
- Freeze source-label inclusion, normalization, training inputs, baseline
  definitions, metric, aggregation, margin and failure handling before exposing
  treated Parse values to model selection. The existing 72,446-row selection
  includes both B Naive and B Intermediate/Memory; do not choose a favorable
  subtype after scoring.
- Establish donor/participant overlap status as far as source provenance allows;
  retain unknown overlap as a limitation of independent validation.

Count verification necessarily reads treated counts. Keep that engineering work
separate from outcome-driven model tuning, and record when outcomes become
available to the predictive analysis.

## What a completed prediction test could support

Compare a frozen candidate with no-change and training-only mean-response
baselines for every donor. Report donor-level errors and failures as well as the
aggregate, and retain the established 5% improvement gate against both baselines
if claiming continuity with the existing transfer experiments. Report uncertainty
only if its calibration is separately evaluated.

A passing result would support bounded prediction of average RNA response under
this cohort and protocol. It would not validate cell-state trajectories, protein
response, immune recognition, disease progression or clinical benefit. The current
three-study context model remains unpromoted because two held-out studies fail
its existing transfer gate.

The frozen per-cell normalization replay and terminal review now pass for all
twelve donors; see the [retained evidence](../../StreamedLogCPM/ParseReplay/README.md).
The next prediction work is experimental provenance resolution and a documented
admission decision for the unpromoted candidate. The versioned
feature contract and renewed development evaluation are complete; their failed
transfer gates remain part of the admission decision.
Do not tune against Parse treated outcomes to resolve these requirements.

The [versioned exact 11,600-feature panel](../../StudyContextKernel/ExactParsePanel/README.md)
now has a frozen metadata-only contract and a complete renewed three-study
development evaluation. HIRISA passes the unchanged gate; Kang and GSE181897
fail. The candidate remains unpromoted, the original model is unchanged, and
no Parse prediction was fitted or scored. Dose/reagent and participant-overlap
provenance remain unresolved.

## Annotation reference lead inspected 2026-09-12

The [Allen Institute downloads](https://apps.allenimmunology.org/aifi/resources/parse-10m-cytokines/downloads/)
provide AIFI L1/L2 CellTypist models trained from labeled Parse PBS controls and
DESeq2 pseudobulk results across twelve subjects. These are concrete candidates
for reference-mapping interoperability and a method comparison. Their presence
is not evidence that NumiVivo has imported, executed or matched them. No model
or outcome file was downloaded in this inspection.

Because these models use this dataset's PBS controls, evaluation on the same
Parse donors would be a compatibility/development check, not independent donor
validation. Preserve the existing literal B-cell selection; do not substitute
new labels after seeing predictive scores. Before execution, record model file
hashes, feature/preprocessing requirements and reference cohort provenance.
The published DESeq2 results use different labels; direct numerical agreement
requires matched cell membership, filtering and design, not just gene names.

Subsequent [CellTypist execution](../../../ReferenceMapping/CellTypist/README.md)
now qualifies imported-model numerical agreement and native H5AD CSR/CSC
orchestration on 267 retained real cells. It does not change the B-cell selection
or supply independent annotation or prediction validation.

The [Parse experimental description](https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/)
confirms 24-hour exposure but does not specify IFNB dose/reagent in the inspected
text. Dose was unresolved at that inspection; the later primary-workbook
receipt below resolves it.

## Prediction endpoint audit: 2026-09-12

The executed [twelve-donor audit](endpoint-audit.json) inspected the current
native report schema and bound every report hash across all 72,446 cells.
The reports contain integer pseudobulk sums and cell QC, but no per-cell
log-normalized feature means. The context model requires the latter:
mean over cells of log1p(1e6 × count / full-source cell total).

Taking log1p after summing counts would change the endpoint. The next concrete
implementation is therefore a bounded count-stream accumulator for the frozen
per-cell transformation, with independent sparse-reference verification and
explicit treatment of zero-count cells. Preserve all 40,352 features in library
totals before projecting to the frozen 11,600-feature panel. This preparation
was the next step at the time of this audit. The subsequent normalization replay
and exact-panel prediction test are complete; the prediction test failed, as recorded
in the current status below. This historical audit itself fitted no prediction.

The read-only audit script records its original host paths and requires a new
output directory. Its receipt checks are an endpoint inventory, not a replacement
for the complete count/source replay verifier.

## Primary donor metadata located

The authors' [donor workbook at a pinned revision](https://github.com/theislab/HumanCytokineDict/blob/6f9bc00381227fe8b1aa8dcb4e3f6f168cfc3229/annotations/donor_metadata.xlsx)
contains twelve donor rows and source donor identifiers. The
[provenance review](donor-provenance-review.json) binds its SHA256 and the
training cohort metadata inspected. Different identifier spellings or numerical
values across studies do not prove that participants differ. A verified
cross-study identifier namespace or supplier provenance is still needed before
claiming participant independence. No demographic or medical attributes were
used for model fitting or matching. Dose/reagent was unresolved at this point in
the investigation and was resolved by the subsequent workbook receipt below.

## Primary dose and reagent resolved

Europe PMC's supplementary-file service supplied the original `media-2.xlsx`
workbook for PMC12724453 after earlier browser retrieval attempts failed.
Sheet `2.cytokine_screen`, row 43, identifies IFNB1 / IFN-beta at **100 ng/mL**,
**Biotechne / R&D systems**, catalogue **8499-IF-010/CF**. The
[dose receipt](dose-provenance.json) records the exact header/row and workbook
SHA256. The complete workbook is retained on the execution host; response-analysis
sheets were not used for model selection. This supersedes earlier missing-dose
statements in the historical audit above.

Dose/reagent identity is now known. Cross-study exposure equivalence and
participant independence remain unverified. Neither this metadata result nor
the completed normalization replay overrides the candidate's failed development
transfer gates. A Parse evaluation must preserve the frozen model, all donors,
full-axis normalization, both simple baselines and the existing improvement
criterion; it cannot be presented as protocol-matched independent validation
without the additional provenance. No Parse prediction has been fitted or scored.
