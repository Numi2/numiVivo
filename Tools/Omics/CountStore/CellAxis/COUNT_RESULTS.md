# Complete same-cohort count qualification, 2026-09-11

**Both full native ingestion and original-source replay pass.** The file-backed and resident owners consume the same complete canonical stream for all **725,031 Parse PBS/IFN-beta cells, 12 donors, 24 groups, 40,352 features and 1,373,870,697 records**. Each phase checks all 3,456 original source runs against the earlier donor qualification; no downsampling or new QC exclusion is applied.

All original identity bytes, cell QC values, group memberships and aggregate coordinates agree between owners, prior retained counts and fresh independent sums. The retained total is **3,070,817,047 counts**. Both phases send 21,981,931,152 bytes with SHA-256 `caaeb0e4ffbf0775764ebf918438f9d3c7c4c536920e8224db5ade5c6b59a6e7`.

| Phase | File peak RSS (MiB) | Resident peak RSS (MiB) | File seconds | Resident seconds |
| --- | ---: | ---: | ---: | ---: |
| ingest | 165.6 | 1005.0 | 2910.1 | 2912.5 |
| verify | 181.8 | 1024.7 | 2834.1 | 2835.3 |

These are per-native-child process measurements. Elapsed time includes waiting for the bounded remote source stream and does not measure isolated compute throughput. Peak memory is lower for this complete count path; it does not establish bounded memory for every downstream algorithm or every admitted maximum. The [initial full statistical analysis](../Expression/README.md) peaked near 1 GiB; its [incremental report writer](../Expression/Memory/README.md) now reduces the same complete run to 202.3 MiB. The statistical model still retains per-feature state.

The count executable is SHA-256 `cae3f6e173b4a9c52f6f606743d99bfe532cad76b69f3cccd83dacb8eee18f9a`, built from 105 frozen inputs at source commit `05c01ba1a8f4fc4efcefa5df7684a51a9f7f7444`. Its 14 cell-axis/count tests pass. The prior executable is `20b13e585e1526dba7be566484373a3e4775efebd66f22b3b20dff31c1ee3516`. These count results remain tied to those binaries; the later expression executable is a separate qualification.

The original historical QC discrepancy remains: source totals exceed retained matrix totals by 32,500 counts and detected-feature totals by 31,902; each QC field differs in 30,634 cells. Matching retained matrix counts does not erase this discrepancy or prove that both mismatch sets contain the same cells. Source feature-name and dose gates still prevent the separate duration-prediction validation. No predictor is fitted or scored here.

## Retained evidence

[Complete results](evidence/2026-09-11-counts/summary.json), [archive identities](evidence/2026-09-11-counts/manifest.json) and [restoration](evidence/2026-09-11-count-restoration/restoration.json) retain both native bundles, every phase's range/QC records, independent matrices, all process metrics and frozen checkers. The complete offline checker runs after original execution and again against restored artifacts. Restoration verifies archived bytes and arithmetic; it does not perform a third native count run or network replay.

`check_counts.py` checks the dated complete evidence without network or subprocess access. `retain_counts.py` rejects incomplete phases before creating a successful result. `restore_counts.py` combines the exact axis and count archives and uses separately verified preparation/count dependencies. Optional APFS reuse requires matching archive hashes and separate inodes. `finish_counts.py` binds an already-running source controller's identity, then retains and restores only after terminal success; it does not restart failed jobs.

See [dependencies.json](evidence/2026-09-11-counts/dependencies.json) for all required source artifacts. The count archive contains selected source derivatives, not a local copy or whole-file hash of the 227 GB upstream H5AD. Parse data and derivatives are **CC BY-NC 4.0**.
