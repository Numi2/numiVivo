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

A complete library/CLI build and Parse source-stream integration remain to be
executed. No Parse prediction is fitted or scored, and this numerical result
does not establish biological predictive accuracy.
