# Native HNSW qualification, 2026-09-09

This archive binds native builds, full-cohort graph construction, exact-query
recall and independent numerical checks to executable SHA256
`122f99972d2018b1741d614d499c9ef842e45093e2f160f78e9d778f5e201161`.
`runtime.json` records the M4 Pro/24 GiB/macOS 26.6 host, HDF5 hash and native
source hashes verified equal between the local publication tree and remote build.
The base is `e1a8206194ad90017473f2bedc4e746e0857bf8c` plus these recorded changes.
Script hashes are in `checks.json`; vendored upstream provenance is in the source.

## Protocol and results

`protocol.json` was declared before execution: M=16, construction width=200,
search width=128, seed=7, k=15 including self. Mean strict recall must be ≥0.95
and fifth-percentile strict recall ≥0.85. Both gates passed without tuning.

| Cohort | Total cells | Exact-query checks | Mean strict recall | Fifth percentile | Minimum |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baron | 8569 | 8569 | 0.9999666572195456 | 1.0 | 12/14 |
| Hagai | 13863 | 13863 | 0.9999896950773385 | 1.0 | 13/14 |
| Norman | 111445 | 2048 | 0.9994070870535714 | 1.0 | 12/14 |

Small-cohort recall uses the full archived, independently qualified exact graphs
from `../2026-09-09-windowed-neighbors/real/{baron,hagai}/workers-4/`.
Current score bytes are equal to those exact baselines. Norman queries are a
uniform sample with seed 20260909; each is searched with SciPy cdist against
every cell, with explicit distance/index cutoff-tie ordering. This is sampled
Norman recall, not exhaustive all-cell recall. `reference/*/recall-queries.npz`
retains row IDs, exact neighbors/distances and per-query strict/tie-aware recall.
Tie-aware recall is reported separately and does not replace the strict gate.

Every returned distance is checked independently with NumPy; maximum error is
3.56e-15 across cohorts. umap-learn computes the fuzzy graph from the returned
approximate neighbors, with elementwise tolerance `1e-5 + 1e-4*abs(reference)`.
Topology matches exactly for all three; maximum weight error is 7.57e-6.
Bandwidths, kernel mass residuals, symmetry, isolation and components also pass.
This fuzzy check does not claim exact-neighbor membership. Baron/Hagai fuzzy-edge
Jaccard against exact graphs is 0.9999434575/0.9999727132.

## Native lifecycle and memory

| Cohort | Fit seconds | Graph publish seconds | Publish peak resident bytes | Verify seconds | Verify peak resident bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baron | 6.67 | 6.85 | 222674944 | 6.89 | 229638144 |
| Hagai | 13.35 | 15.11 | 300728320 | 15.18 | 301973504 |
| Norman | 135.86 | 161.56 | 1357234176 | 162.54 | 1353629696 |

These are single-run native CPU lifecycle measurements including PCA
reconstruction, HNSW, fuzzy graph and serialization. External reference checks ran
on the laptop; no same-host scverse speed comparison or Metal claim follows.
`real/commands.json` retains all nine successful timed commands and their logs.
Graph verification reconstructs both input and graph with the current executable.
Norman performs 467109415 metric evaluations including repeats; its decoded cache
peaks at 17831200 bytes and estimated serialized index size is 16998788 bytes.
Cache/index counters exclude allocator and container overhead. Graph JSON and
arrays remain resident and explain why overall RSS exceeds cache/index counts.
The three real cohorts fit in the default score cache. The separate 600×64 Swift
regression forces eviction with a two-tile cache and checks exact results/replay.

All 75 Swift tests in 22 suites passed. The full release and scoped H5AD build
passed. The first scoped invocation failed with exit 126 because its script is
not executable; its log is retained. Reinvoking through `bash` succeeded. There
was no recall-gate failure. `fixtures/` retains 20 HNSW commands, including ten
expected rejections, fitted/query input coverage, deterministic replay, smaller
cache invariance, small-fixture exact results and rehashed tampering rejection.
`fixture-inputs/` retains 22 query-input regression commands with 13 expected
rejections. A fresh full-Baron exact-mode run also matches prior graph SHA256
`b2d753af1463ed6f07b7b892bb52b592a762848787f23860602439fcf5889df1`.

## Reproduction and scope

Decompress `.gz` files before replay. `real/*/graph/` preserves the graph, execution
report, receipt and full PCA state except the original H5AD. Source H5AD hashes
remain in input receipts; source paths and fit plans are in the command manifest.
Original data provenance is documented in the existing Baron/Hagai benchmarks
and full-Norman PCA evidence. Native replay requires those original sources and
the qualified runtime. A different binary must refit from source, not rehash old
receipts. `inputs.json` records the original runner inputs; archived input plans
are under `real/CASE/graph/input/plan.json` if those original paths are unavailable.

Run `run_hnsw.py --binary BINARY --inputs INPUTS --out NEW_DIRECTORY` with
`NUMIVIVO_HDF5_LIBRARY` set to the qualified HDF5 library. Run
`check_hnsw_reference.py --bundle GRAPH --case CASE --protocol protocol.json
--out NEW_REFERENCE`, adding `--exact-graph` and `--exact-scores` from the exact
baseline for Baron/Hagai. The checker accepts gzip baseline files. Its Python
versions are recorded with each result. The global `SHA256SUMS` binds archive bytes.

This milestone qualifies the declared approximate search and numerical graph
construction. It does not qualify million-cell graph storage, downstream
clustering/embedding biology, donor integration, unseen perturbations, protein
or reaction kinetics, tissue outcomes, or an end-to-end Metal speed advantage.
