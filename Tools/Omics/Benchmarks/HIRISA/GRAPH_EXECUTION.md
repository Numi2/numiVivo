# Frozen full HIRISA neighborhood graph qualification

PCA publication `e4b630b6` qualified all 1,612,594 cells, 18,082 genes and
32,251,880 scores. This next experiment uses those same 20 PCs and every cell;
it does not refit preprocessing on a held-out split or qualify predictive biology.

## Resource admission before execution

The C++ bridge still rejected more than one million rows after Swift PCA/storage
admission was raised to two million. Its independent 128-million-neighbor-entry
cap also disagreed with the supported row and k bounds. Shared C/Swift constants
now own the row/dimension limits; streamed output admits the corresponding
product. The four-million-entry resident JSON cap is unchanged. Graph snapshot
and verification byte bounds now derive from the maximum symmetric graph size.

The previously qualified 111,445-cell Norman run used 467,109,415 metric calls.
Linear extrapolation alone exceeds the former two-billion supported ceiling on
HIRISA. The new supported ceiling is 32 billion calls, still checked on every
metric evaluation. The default remains 500 million. Cache support extends from
64 MiB to 1 GiB, while the default remains 32 MiB. This permits the complete
258,015,040-byte 20-PC representation to fit without repeated tile eviction;
it does not make the HNSW index out of core.

A focused native test creates a sparse score file for 1,000,001 rows with k=128
and reaches the work-budget rejection through the real C++ index. It also checks
the two-million-row boundary and oversized cache/work rejection. This is resource
contract evidence, not a synthetic substitute for the full experiment.

## Frozen experiment

- Complete GSE306664 HIRISA source; source order and identities unchanged.
- Scores SHA256: `62e6b0644e2c954749b77946ba748e184c33580614f5756b1bc1b0980df86296`.
- Binary graph output; k=15 including self; local connectivity=1.
- HNSW connections=16, construction width=200, search width=128, seed=7,
  serial insertion/query. These are the prior qualified search settings.
- Explicit maximum distance evaluations=16,000,000,000, more than twice the
  linear Norman work extrapolation. This is a resource allowance, not a claim
  that work will be linear or that the graph will fit the allowance.
- Explicit score cache=258,048,000 bytes, rounded up to complete 256-row tiles.
- Uniform fixed sample of 2,048 queries, NumPy `default_rng(7)`, sorted in source
  order. Each query is searched exactly against **all 1,612,594 candidates**.
- Exact FP64 SciPy Euclidean distances; cutoff ties resolved by source row.
- Mean strict recall ≥0.95 and fifth-percentile strict recall ≥0.85. Tie-aware
  recall is reported separately and cannot replace either gate.
- Every returned distance checked against NumPy, rtol/atol=1e-12.
- Every bandwidth and fuzzy edge checked against umap-learn on the returned
  approximate neighbors, with the existing rtol=1e-4, atol=1e-5 weight gates.
- Every binary coordinate, source score/identity/QC hash, graph symmetry,
  connected-component count and isolated-cell count checked independently.

Frozen protocol SHA256:
`f48c1a42a16608d80d372f999ab77c9e6cd20d1a03d07409d1b9a0b15dcda6fa`.
The protocol and plans are retained beside the execution receipts. The exact
query panel is generated before native graph output is inspected.

## Execution and evidence order

1. Focused native HNSW/PCA/binary-store/window tests, then complete release build.
2. Freeze executable, authored source hashes, HDF5, compiler and hardware.
3. Fitted/query graph lifecycle, tampering and resource controls; full Baron and
   Hagai graph regression against previously independently qualified records.
4. Refit complete HIRISA PCA using that executable (older executable receipts
   cannot be migrated by rehashing). Require all six fitted artifacts to equal
   the already-qualified PCA bytes.
5. Publish complete binary graph; reconstruct it with the same executable.
6. Run `reference_graph.py check` against the frozen exact panel and full graph.
   Retain failures before changing settings. Archive complete evidence and
   publish only once native replay and independent checks are both established.

`reference_graph.py prepare` holds the compact PCA matrix and at most sixteen
query-by-cell distance rows. `check` uses mapped native records and chunked
returned-distance checks; the independent scverse fuzzy graph remains a resident
sparse reference. Native and reference hosts differ, so these timings do not
establish a controlled native/scverse speed comparison. Graph recall and numerical
agreement do not establish clustering, batch integration, biological signal
preservation, or Metal acceleration. Those remain subsequent gates.
