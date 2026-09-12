# Streaming per-cell log-CPM means

`VivoStreamedLogCPM` computes group means of log1p(1e6 × count / full-cell total)
from ordered positive count records. It retains group-by-feature sums and
compensation terms, plus cell metadata; it never constructs a cells-by-features
matrix. Full-axis cell totals are required and verified against incoming records
before a result can finish. Projection belongs after normalization.

Zero-count cells contribute zero and remain in the mean denominator. Duplicate,
descending, out-of-axis, zero-valued, excessive and incomplete records fail.
An error closes the accumulator; neither further records nor a partial result
can then be accepted. Every declared group must contain cells.

The API limits groups to 4,096, cells to 10 million and group-by-feature entries
to 16 million. These admission bounds are not million-cell performance evidence.
It is included in the scoped library source list. There is no new CLI or Parse
orchestration yet, and existing count-stream receipts remain unchanged.

## Executed verification

The actual Swift source compiled with Swift 6 optimization on the Mac mini.
All 2,711 cells, 36,601 RNA features, 5,218,473 positive records and 11,786,194
counts from the existing 10x source were consumed. Three deterministic row
partitions exercise grouped computation; they are not biological annotations.
Independent SciPy sparse log-normalization and means agree within 1.96e-14,
against a frozen numerical tolerance of 1e-11. Zero-cell, all-zero, closed-state
and malformed-record checks pass.

The first test harness compile failed because a throwing call was inside a
nonthrowing precondition autoclosure. Moving that call outside the assertion
fixed the harness; the production source required no change.

The report binds the source, canonical stream, plans, reference and native
outputs. Scripts retain their original remote paths. Adjust paths in a fresh
study to reproduce; prepare.py refuses to overwrite counts.bin and the Swift
harness refuses to overwrite native.json. Build the harness with the actual
Sources/NumiVivoKit/Omics/VivoStreamedLogCPM.swift file and Main.swift using
`swiftc -swift-version 6 -O -parse-as-library`, then run prepare.py, the harness
with the study path, and verify.py. NumPy, SciPy and h5py are required.
The retained remote study is /Users/n/numivivo-streamed-logcpm-20260912.

The complete scoped library/CLI build at 1a2f6bba passed. Linking the harness
against that library reproduced every standalone result value exactly and
passed the malformed/zero-cell checks again. product-verification.json binds
the library, CLI, source inventory, stream harness and product output. Existing
Metal API deprecation warnings remain in the build log.

The [twelve-donor Parse normalized-mean replay](ParseReplay/README.md) has completed
in /Users/n/numivivo-parse-logcpm-20260912 and passed terminal review.
run_parse.py uses the previously qualified adapter, verifies source range and
complete canonical-stream hashes, retains all frozen B cells, and compares every
group/feature mean to independent NumPy accumulation at absolute tolerance 1e-10.
Its protocol binds the driver and runtime before execution. StreamMain.swift
uses the actual product library and handles records split across pipe reads.
The shared source controller lock prevents overlap with an exclusive ingestion.

No Parse prediction is fitted or scored, and these numerical checks do not
establish biological predictive accuracy. Full-app build and end-to-end performance remain outside this result.

## Terminal review

verify_parse_complete.py independently inventories all twelve donor plans,
cell denominators and zero counts, 40,352-feature axes, retained source range
receipts, canonical stream identities, native outputs and reference arrays.
Supply --driver-sha256 from the trusted published run_parse.py, not a digest
chosen by the study itself. It rejects incomplete cohorts and mismatched inputs.
It rechecks retained evidence; it does not independently refetch source counts.

The live incomplete study was rejected, and four admission checks passed:
missing completion, forged completion without donors, changed tolerance, and
driver mismatch. The positive complete-cohort review now passes for all 72,446 cells and
40,352 features, with maximum absolute error 7.99e-15. test_parse_complete.py uses a separate temporary directory
and does not change the running study.

The [native CLI](CLI/README.md) now exposes the shared bounded count-stream
consumer with required SHA256 verification and output protection. Real-matrix
and malformed-input checks pass; this adds no biological prediction claim.

The [feature-stream bundle](FeatureStream/README.md) emits one feature at a time
and qualifies the full Norman result at 44.1 MB peak CLI RSS.
