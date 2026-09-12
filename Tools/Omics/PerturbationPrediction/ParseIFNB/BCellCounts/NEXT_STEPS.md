# From validated B-cell counts to a predictive test

The running count workflow prepares a possible external RNA-response test.
Completion establishes faithful count processing; it does not establish predictive
accuracy. Leave its frozen inputs and live driver unchanged while it runs.

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

- Resolve the original IFN-beta dose and reagent from an accessible primary
  source. The established 24-hour duration does not resolve dose comparability.
- Resolve the feature contract: the current 11,800-feature context model has
  11,600 exact source matches and 200 absent features. Do not silently zero-fill,
  alias or choose a panel after inspecting outcomes. A changed panel requires an
  explicitly versioned model and renewed development validation.
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

While ingestion runs, useful work is limited to checks, provenance resolution and
protocol preparation that do not alter the frozen execution or tune against its
treated outcomes. Another heavy fit is not a prerequisite for these steps.
