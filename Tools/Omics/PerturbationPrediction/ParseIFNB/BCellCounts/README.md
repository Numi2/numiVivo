# Complete Parse B-cell count workflow

See [the final evidence review and prediction admission steps](NEXT_STEPS.md) for
the work that can proceed during ingestion and the remaining scientific gates.

The frozen cohort contains **72,446 source-labeled B cells, 124,909,573 stored
count records and all 40,352 RNA features** across twelve PBS/IFN-beta donor
pairs. It is the exact metadata selection in [BCellAdmission](../BCellAdmission/README.md).
No donor, selected row or feature is removed on QC disagreement or response.

All twelve native ingestions and source replays passed. The terminal review
validated every donor bundle, independent cell and donor/condition counts,
source range hashes and canonical stream hashes. The complete cohort contains
**253,128,870 total counts**. No prediction was fitted or scored.

The [complete evidence archive](numivivo-bcell-complete-20260912.tar.gz) and
[member manifest](numivivo-bcell-complete-20260912.tar.gz.manifest.json) retain
13,990 files, including the exact native binary and source dependencies.
The 33,834,827-byte archive passed independent member verification remotely and
again after transfer. SHA256:
`5309bf51927b7abec526257cac948af899b67996c47be6d1859a24590a4d9e37`.
This qualifies selected count handling, not biological predictive accuracy.

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

Run the read-only terminal review with the repository's preparation manifest:

```sh
python verify_complete.py /path/to/study --preparation-manifest manifest.json
python -m unittest test_verify_complete.py
```

Keep the verifier outside the running study so its frozen driver files remain
unchanged. It exits nonzero with `not-verified` for incomplete or inconsistent
evidence. The admission tests cover missing completion, altered preparation,
wrong cohort totals and escaping paths. The live incomplete remote run was also
rejected. The full successful-cohort terminal review now passes; biological prediction
remains unqualified.
The verifier reuses the retained independent bundle checker and rejects Python
optimization mode because that checker requires assertions. Use a trusted
published preparation manifest; a self-authored manifest is not a provenance
anchor. Archival packaging and independent archive verification remain separate.

After terminal review, `package_complete.py` takes a shared controller lock,
reruns the review and streams every study file into a compressed archive outside
the study. Only the controller lock and Python bytecode caches are excluded;
failed attempts, source dependencies and the exact native binary are retained.
It checks the source inventory again and independently reads every archived
member to verify exact membership, sizes and hashes before finalizing the archive.
Failed packaging retains its `.partial` file and never overwrites an existing
destination. A sidecar records the reviewed inventory and archive hash.

```sh
python -m unittest test_verify_complete.py test_package_complete.py
python package_complete.py /path/to/study --preparation-manifest manifest.json --output /separate/path/complete.tar.gz
```

Six software tests passed. A real invocation also rejected the active remote
controller without creating an archive. Complete-cohort packaging and independent member verification passed; archive
integrity alone does not establish biological accuracy. Keep both Python tools together outside the study.

verify_archive.py checks this preparation archive. It embeds the twelve plans,
their row maps, scripts and software verification. Original source axes/chunk
maps and the native binary remain in the hash-bound published dependencies,
rather than being duplicated here. It does not embed raw count streams or
completed biological validation. The retained scripts' original absolute paths
describe their execution host and must be adjusted on another host.

The IFN-beta dose/reagent and model feature contract remain unresolved.
No prediction is fitted or scored by this workflow. Original Parse metadata and
count derivatives retain **Parse Biosciences, CC BY-NC 4.0** attribution.

## First completed donor: partial execution evidence

The later [complete ingestion review](ingestion-review.json), executed by
[review_ingestion.py](review_ingestion.py), checks all twelve native bundles,
their planned identities and independent counts, per-run QC sums, stream/receipt
hashes and the complete ingestion receipt. It covers 72,446 cells, 124,909,573
records and **253,128,870 total counts**, with 3,433 source runs and 6,866 range
receipts. Historical source gene-count and transcript-count QC each disagree
for 3,188 cells; their aggregate excesses are 3,305 detected features and 3,345
counts. No affected row was removed. The source and preparation hashes match.

This historical ingestion review preceded replay. The complete archive above
now embeds every donor bundle and replay receipt. No prediction was fitted or scored.

[The retained offline review](donor1-ingest-review.json) rechecked Donor1's
native bundle against independent counts: 4,649 cells, 9,621,618 records and
22,464,566 total counts across 288 source runs with 576 range receipts.
All preparation members and bound dependencies matched their published hashes.
Historical source gene-count and transcript-count QC each disagree for 293
cells; their totals exceed the selected matrix by 306 detected features and
310 counts, respectively. These discrepancies were retained, with no row removal.

[The review script](review_donor1.py) records its original remote paths and
requires the retained study dependencies. This is an offline ingestion review
of one real donor, not a second source replay, a complete-cohort result, or a
successful run of the full terminal verifier. The later complete terminal review above supersedes this partial status.
