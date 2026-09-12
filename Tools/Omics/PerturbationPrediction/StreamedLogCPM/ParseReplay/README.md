# Complete Parse normalized-expression replay

All twelve frozen donors pass native streaming accumulation and the terminal
review: 72,446 B cells, 40,352 features, and both IFN-beta/PBS conditions.
Maximum absolute error against independent NumPy per-cell accumulation is
7.993605777301127e-15, below the pre-execution tolerance of 1e-10.
No prediction was fitted or scored.

The endpoint is mean per-cell log1p(CPM), using each cell's full source feature
count total and retaining zero cells in group denominators. It is not a
log-normalized pseudobulk sum. Source range hashes and canonical count streams
match the previously verified count ingestion for every donor.

[terminal-review.json](terminal-review.json) records the complete-cohort review.
[evidence.tar.gz](evidence.tar.gz) contains all 3,498 retained study files,
including the frozen protocol, driver, donor plans, source range receipts,
native results, independent reference arrays and completion receipts.
[archive-manifest.json](archive-manifest.json) binds every member and the archive.
All member hashes matched before/after packaging and after transfer.

The archive is 14,837,916 bytes, SHA256
`beb2b132c942e6c97b70541344cec65baf03c30875f75412f4729712931033af`.
Original count-source evidence and the qualified product runtime are separate
[dependencies](../README.md); this is not a standalone raw-source reproduction.
The terminal review rechecks retained receipts and arrays, without independently
refetching raw counts. The replay itself refetched and verified source ranges.

This closes the normalization input gap. Dose/reagent and participant-overlap
provenance remain unresolved, and the versioned feature-panel candidate failed
transfer gates in two of three development studies. Numerical agreement does
not establish biological predictive accuracy or end-to-end performance.
