# Binary graph-store qualification, 2026-09-09

Complete Baron (8,569 cells), Hagai (13,863 cells) and Norman (111,445 cells)
were fitted with 20 PCs and processed with k=15 HNSW using the previously declared
M=16/construction=200/search=128/seed=7 configuration. Complete Baron also exercises
four-worker exact search with a 200-million unique-pair budget. No cells are
subsampled for binary graph comparisons.

| Cohort | Binary publish seconds | Binary peak resident bytes | JSON publish seconds | JSON peak resident bytes | Binary verify seconds | Verify peak resident bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Baron | 6.82 | 202817536 | 6.83 | 229736448 | 6.85 | 204668928 |
| Hagai | 14.80 | 226689024 | 14.88 | 310018048 | 14.87 | 224559104 |
| Norman | 157.97 | 477413376 | 161.43 | 1420509184 | 160.15 | 468582400 |

Norman's observed resident peak is 66.4% lower with binary storage. Full-Baron
exact binary publication/verification takes 6.34/6.35 seconds with resident peaks
210468864/197607424 bytes. These are single-run native CPU lifecycle measurements
on Apple M4 Pro/24 GiB/macOS 26.6, including input PCA reconstruction. Completed
files were copied for local independent checks during later cases. There was no
competing native build or benchmark. These timings do not establish general speed
improvements, GPU execution or a same-host scverse comparison.

Every neighbor ID/distance, rho, sigma, mass residual, CSR offset, edge coordinate
and FP64 weight bit matches the previously independently qualified full graph.
All PCA score bytes and cell identities also match. Baron exact uses the exact
reference in `../2026-09-09-windowed-neighbors`; HNSW uses `../2026-09-09-hnsw`.
`reference/*/checks.json` binds each exact baseline by graph/score SHA256. The
same-executable JSON graph payloads are byte-identical to those HNSW baselines,
as recorded in `json-baselines.json`.

Norman has 1,671,675 fixed-width neighbor entries, 2,513,474 connectivity entries
and one connected component. Its four binary files total 74,094,880 bytes. The
binary equality check covers every record. Recall retains the earlier scope:
99.9407% mean strict recall on 2,048 preselected queries searched exactly against
all 111,445 cells. Equality does not turn sampled recall into exhaustive recall.
Baron and Hagai retain all-cell mean strict recall 99.9967% and 99.9990%.

## Runtime and publication boundaries

The full-cohort benchmark executable SHA256 is
`bd93e9330383baf1666a2c77ce7d2e4a53f8d7e5ed51b40e7c46920332f9bd6a`,
built from `85e2f23a0810bddd2a86b6c68a3b241d2683a3b0` plus the source changes
fingerprinted in `runtime.json`. Local publication and remote build source bytes
were compared before and after validation. HDF5 2.2.0 is bound by its exact hash.
All 77 Swift tests in 23 suites passed (66 single-cell, nine multi-assay and two
count-store tests). The full release build passed in 228.21 seconds. The scoped
build produced a fresh executable that passed HDF5 import/export, with all fixture
counts and original source bytes independently checked.

During the benchmark, separate HCC1395 workflows and genomic CLI routing changes
advanced upstream to `b359e71182bda52cf13fbfaddd22ca71b22b5a46`. Those changes were
preserved with a fast-forward merge. The only upstream compiled-source change was
`Sources/NumiVivoCLI/main.swift`; its new handler selects explicit genomic commands
and does not claim single-cell commands. All omics/C++ source bytes remained equal
to the benchmark's hashes. The merged release build passed in 11.72 seconds and
has SHA256 `e3fb4fc40ee336a5c7b7d1e538d52f7e10788e5b52cf9dedc36f503d9e0ded60`.
Fresh PCA/query and binary-store CLI suites passed with this merged executable.
`publication-scope.json` records both identities. The large-cohort timings and
receipts remain attributed to the benchmark binary; they were not relabeled as
runs of the merged executable.

An intentional mixed-binary input attempt was rejected with exit 65. This is an
expected provenance control, not an unresolved failure. There were no failed
native benchmark phases or numerical comparisons in this milestone. Expected
negative-control outputs and existing compiler warnings remain in the archive.

## Artifacts and replay

`real/commands.json` retains all 14 native command arguments, statuses, wall times
and resident peaks; neighboring logs preserve output. `real/CASE/binary/` retains
the graph metadata, four binary files, plan, receipt, execution report and PCA
input state. `real/baron/exact-binary/` retains the complete exact comparison.
JSON graph/plan/receipt/execution outputs are also archived. Where JSON input is
omitted, its input receipt equals the sibling binary input receipt; use that same
input bundle for reconstruction. Real original H5AD files are omitted and remain
bound by input receipt hashes and source paths in the command manifest. Dataset
provenance is in the existing Baron/Hagai and Norman benchmark documentation.

`fixtures/` and `merged-fixtures/` each contain 22 graph commands with six expected
rejections: fitted/query inputs, exact/HNSW searches, every binary field against
JSON, repeatability, overwrite/format rejection and rehashed changes to all four
files. `fixture-inputs/` and `merged-fixture-inputs/` each contain 22 PCA/query
regression commands with 13 expected rejections. The latter uses fresh inputs from
the merged binary, not earlier executable-bound receipts. `scoped-smoke/` is a
16-cell numerical fixture and is not biological evidence.

Decompress `.gz` files before replay. `run_graph_store.py` accepts a binary, an
input manifest with name/source/fitPlan fields and a new output directory. Sources
and plans can be recovered from the native command manifest and archived PCA input
plans. Set `NUMIVIVO_HDF5_LIBRARY` to the fingerprinted runtime.
`check_graph_store_reference.py --bundle BUNDLE --baseline-graph GRAPH.json.gz
--baseline-scores scores.bin.gz --out NEW` checks every binary record using NumPy
and SciPy. Its baseline graph comes from the earlier exact/HNSW qualification;
those references already include independent distance/fuzzy-graph checks. Archive
checks and scripts are fingerprinted in `checks.json` and `SHA256SUMS`.

The graph is now stored and constructed through bounded file windows and a disk
transpose/merge. Input PCA state, identities, HNSW index and cell-scale degree,
offset and component arrays remain resident. Native file-backed clustering,
embedding and integration, actual million-cell qualification, downstream biological
stability and Metal/scverse performance comparisons remain open.
