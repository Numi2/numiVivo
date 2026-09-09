# Windowed exact-neighbor qualification, 2026-09-09

Complete Baron (8,569 cells, 20 PCs), complete Hagai source-QC cohort (13,863
cells, 20 PCs), and the complete held-out Baron human3 query (3,605 cells,
20 PCs) were processed. Each uses k=15 including self and an explicit 200-million
unique-pair budget. No cells were subsampled for these graph comparisons.

| Input | One-worker publish seconds | Four-worker publish seconds | Four-worker maximum resident bytes | Verify seconds |
| --- | ---: | ---: | ---: | ---: |
| Baron | 7.15 | 6.61 | 232783872 | 6.68 |
| Hagai | 15.92 | 14.66 | 303644672 | 14.68 |
| Held-out human3 | 5.08 | 5.02 | 170164224 | 5.01 |

These are single-run end-to-end native CPU measurements on the Mac mini; they
include input/PCA reconstruction, graph construction and output serialization.
Both worker settings perform the same directed-distance work. They are not
isolated kernel measurements or a same-host scverse performance comparison.
Source PCA reconstruction dominates the elapsed time and memory in these runs.

All one/four-worker graph bytes are identical. Independent SciPy checks every
neighbor identity/order and distance (maximum distance error zero for all three).
umap-learn matches the complete fuzzy topology: 176,850, 293,178 and 74,824
connectivity entries respectively. Maximum weight error is below 5.7e-6 against
its float32 reference. Bandwidths, kernel mass residuals, components and isolation
counts pass the existing tolerances. Input Baron/Hagai PCA numerical fields remain
exact versus previously published results; rebuilt query-score bytes match the
previous frozen-query release.

`real/` retains both worker graph/plan/receipt/execution artifacts, input PCA state
and all attempted native command logs. Source H5AD files remain at the immutable
source paths in commands and `inputs/qualified-inputs.json`, bound by the input
receipts. Decompress `.gz` files before replay. Source provenance is in the
existing Baron and Hagai benchmark documentation and full-cohort PCA evidence.
For the query case, training/query source plans are also archived under `inputs/`.

Run `run_window_neighbors.py` with `inputs/qualified-inputs.json` and a new output
directory to reproduce all three native cases using the qualified executable and
HDF5 runtime in `host.txt`. Run `window_neighbors_reference.py` on each four-worker
bundle; it reconstructs the derived reference JSON, which is omitted here to avoid
duplicating graph and score payloads. `reference/` retains independent checks,
versions and outputs; source/script/binary hashes are in `checks.json`.

Failure history is retained. The initial test build exceeded Swift's type-checking
time for a new fixture initializer; it was split into simple expressions and all
73 tests in 21 suites subsequently passed. The first fixture and real query-input
attempts used an older executable-bound PCA receipt and were rejected. No receipt
was relabeled. Inputs were rebuilt from their retained sources with the current
executable; only the failed query branch was resumed, preserving the already
completed Baron/Hagai runs. `real/real-commands.json` therefore contains eight
successful commands and the stale query failure; `real/real-query-current-commands.json`
contains five successful commands. The final runner adds explicit query-input
reconstruction for single-invocation reproduction.

`fixtures-current/` retains 19 graph commands including nine expected rejections,
fitted/query support, one/two/four-worker and tile-size invariance, repeats and
rehashed graph/work/input tampering. `current-fixture-inputs/` retains the 22-command
query-input reconstruction regression with 13 expected rejections. Both fixture
input kinds also pass the independent graph checker. `initial-stale-fixture/`,
`initial-typecheck.log.gz` and all original real command logs preserve the failures.

Scores use bounded tiles and 16 MiB mappings per worker. Final kNN/fuzzy graphs,
identities and input reconstruction remain resident. Exact search is still
quadratic, with at most 500 million unique pairs; this does not qualify full
Norman/million-cell graphs, approximate recall, Metal, embeddings, clustering,
integration or biological prediction.
