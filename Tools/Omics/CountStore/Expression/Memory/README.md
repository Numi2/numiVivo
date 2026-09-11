# Incremental expression-report encoding

The file-backed expression commands now encode feature records incrementally and hash the same byte stream while writing or verifying it. They preserve the existing JSON report, all model calculations, every diagnostic and the source/analysis implementation boundary. The change removes the whole-report `JSONEncoder` tree and `Data` buffer from this product path.

## Why this owner changed

A native phase probe linked to the prior frozen library recomputed the complete Parse PBS/IFN-beta result. Peak RSS was **187.2 MiB after fitting** and **904.0 MiB after whole-report encoding**. The 108,052,499-byte report's SHA-256 remained `4919f0c96dbe114707dd58b5939d708b6a3e9e7a712cb91380edd2424d92ed11`. This probe isolates phases; the earlier product command, including publication, peaked at 970.4 MiB. It is not a new biological experiment or a fresh simultaneous product timing control.

`VivoExpressionReportJSON` streams the main feature table, NB feature diagnostics and large QL fit/test arrays, retaining canonical Foundation encoding for each record and smaller metadata. Fixed ASCII envelope keys follow the existing sorted schema. Optional fields stay absent or present exactly as in synthesized Codable. A 64 KiB output buffer, incremental SHA-256, cancellation checks and the existing total report-byte limit protect publication and verification. The private staging directory is published only after the report, plan, owned source snapshot and receipt complete.

The statistical model still retains per-feature fits and diagnostics, selected design metadata and cohort summaries. This change does not make every analysis allocation independent of feature count, prove bounded memory at every admitted maximum, accelerate model fitting or qualify Metal execution.

## Qualification

The scoped build includes 108 explicit native inputs. All **37 tests across six suites pass**. Exact-byte tests cover log-linear, NB Wald, likelihood-ratio and adjusted QL reports, resident/file membership, optional fixed/empirical effect priors, retained QL failure diagnostics, Unicode/escaped strings, byte limits, sink errors and cancellation. Existing source-identity, altered-report rejection and native reconstruction tests remain included.

The full-cohort product protocol is fixed before execution: all 725,031 cells, 40,352 features and 12 donor pairs; every byte must equal the previously qualified report; native publish peak RSS must be at most half the prior 970.4 MiB; elapsed time must be at most 1.2 times the previous 376.2 seconds. No genes, cells, model parameters or diagnostics are removed. The six earlier edgeR/limma/DESeq2 comparisons are not rerun because an exact full report match preserves their statistical inputs and results.

Both complete native commands pass, and direct comparison confirms all **108,052,499 report bytes** are identical. [Retained process records](evidence/2026-09-11/pipeline.json) and the [frozen protocol](evidence/2026-09-11/protocol.json) distinguish numerical, memory and elapsed-time gates.

| Native command | Peak RSS (MiB) | Elapsed seconds | Result |
| --- | ---: | ---: | --- |
| Prior product publish | 970.4 | 376.2 | Retained baseline |
| Incremental publish | 202.3 | 365.5 | PASS: exact bytes, memory and elapsed gates |
| Incremental numerical verification | 202.1 | 365.4 | PASS: full reconstruction |

Peak publish RSS is **79.16% lower**. The elapsed guard passes; the small timing difference is not an isolated or replicated speedup claim. The new executable is SHA-256 `5f036fdb4b4b13855d910afec37116de9eab028ed19ddb5e1bf765418349b03f`. Model fitting, all 33,899 tested features, 5,946 support exclusions and 507 low-count exclusions remain unchanged.

## Evidence storage

`run_parse.py` measures actual product publish and reconstruction as separate native children, then compares every report byte directly against the exact earlier archive read from its committed Git blob. This avoids allocating another archive working copy. Numerical equivalence, memory and elapsed-time gates are reported separately.

`retain.py` stores changed objects, executable/source, protocols, tests, phase profile and actual process records. Identical report/source objects reference the exact earlier archive, identified by commit, path and SHA-256. Both archives are required to restore. `retain.py restore --archive ... --repo repository --out fresh-directory` reads that pinned baseline Git object; optional `--study completed-study` permits verified separate-inode APFS copies. All restored bytes are hash-checked. Restoration itself does not repeat native inference or network extraction.

Parse data and derivatives are **CC BY-NC 4.0**. Storage and numerical evidence do not add a prediction success. The [biological assessment](../../../../../Documentation/BiologicalPrediction.md) and unresolved independent prediction inputs remain applicable.

The [new archive manifest](evidence/2026-09-11/manifest.json) binds a 3,598,609-byte archive (SHA-256 `ea7eeac4f67d4949c14f294eb87314b77529cd149514d4b5afb978c1ef1f1ffe`). Its exact prior-archive dependency is declared in `contents.json` and the inner `dependencies.json`; the report and source objects are recoverable even when their working copies have been removed.

[Restoration passes for all 148 members](evidence/2026-09-11-restoration/restoration.json): 32 source APFS copies, four copies within the new destination and 112 extractions, all hashes exact. The two measured product executions above establish numerical behavior; restoration checks the retained bytes and does not claim a third inference run.
