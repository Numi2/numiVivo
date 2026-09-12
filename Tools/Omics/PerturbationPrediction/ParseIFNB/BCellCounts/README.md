# Complete Parse B-cell count workflow

The frozen cohort contains **72,446 source-labeled B cells, 124,909,573 stored
count records and all 40,352 RNA features** across twelve PBS/IFN-beta donor
pairs. It is the exact metadata selection in [BCellAdmission](../BCellAdmission/README.md).
No donor, selected row or feature is removed on QC disagreement or response.

This publication records preparation and the execution recipe. Native ingestion
and source replay were launched; **complete count validation is not yet claimed**.
The metadata-only archive does not become count evidence until every native
donor ingestion and replay finishes with independently matching counts.

## Execution and ownership

The same native count-stream owner qualified on the original complete Parse
cohort consumes canonical records. Its exact retained binary SHA256 is
`20b13e585e1526dba7be566484373a3e4775efebd66f22b3b20dff31c1ee3516`.
This uses an existing scoped production binary, not a new full-app build.

The source reader fetches the existing immutable parent runs with the original
ETag and range checks. It validates raw source counts and columns and sorts each
row without merging duplicate columns. A separately tested selection step keeps
only the frozen B-cell offsets and maps them into each native donor stream.
It preserves selected zero-count rows, all feature coordinates and exact counts.
Source runs may include nonselected cells; their transport and QC records remain
explicit, and only selected records enter the new native consumer.

Four fetch workers and eight prefetched source runs bound queued source data.
Each native donor bundle must match every planned cell identity/cardinality,
independently summed cell total, group membership and donor/condition feature
aggregate. All source QC discrepancies remain in the per-run reports.
A second full source pass must match range hashes and canonical stream SHA256
and pass native replay. Completed donor receipts are revalidated before resume;
failed attempts remain separate. A controller lock prevents duplicate drivers.

## Preparation and software checks

All 72,446 selected parent rows occur exactly once across the twelve plans.
Every cell metadata record, feature axis and expected row cardinality equals
its parent projection. Source axes/chunk maps and the parent plan are checked
against the existing published preparation manifest. The selected row archive
is checked against its published B-cell metadata manifest.

The transport tests verify exact row/column/count remapping, a selected zero
row, and duplicate, descending, out-of-range and size-mismatched selection
rejection. Those are software checks, not biological outcome evidence.

## Reproduce

Restore the existing [Parse preparation dependencies](../COUNT_RESULTS.md)
into a separate parent study, and the [B-cell admission](../BCellAdmission/README.md).
With NumPy and the exact native runtime available:

1. Adjust the study/repository paths in prepare.py and verify_preparation.py;
   run them in a new output directory.
2. Set NUMIVIVO_PARSE_STUDY to that directory and run test_selection.py.
3. Run run_bcells.py. It performs ingest and then verify, retaining per-donor
   attempts and validated restart receipts.
4. Require ingest-complete.json, verify-complete.json and complete.json plus
   independent bundle checks before claiming completion.

verify_archive.py checks this preparation archive. It embeds the twelve plans,
their row maps, scripts and software verification. Original source axes/chunk
maps and the native binary remain in the hash-bound published dependencies,
rather than being duplicated here. It does not embed raw count streams or
completed biological validation. The retained scripts' original absolute paths
describe their execution host and must be adjusted on another host.

The IFN-beta dose/reagent and model feature contract remain unresolved.
No prediction is fitted or scored by this workflow. Original Parse metadata and
count derivatives retain **Parse Biosciences, CC BY-NC 4.0** attribution.
