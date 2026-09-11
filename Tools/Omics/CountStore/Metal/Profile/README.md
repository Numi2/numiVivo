# Profile-guided sparse record writer

The shared count writer now uses one initialized, privately owned 1 MiB buffer instead of per-record Swift Array appends. Hashing reads that buffer directly; synchronous file writes borrow its bytes. Raw counts, numerical transforms, record order, admission limits and receipt formats are unchanged. CPU remains the default. This improves repeated normalization timings on the measured cohort, with a slower first CPU call retained below.

## Why this change

Sixteen sampled complete original Kang normalizations matched the retained receipts. Source reconstruction accounts for 52.8% of CPU and 57.3% of Metal main-thread normalization samples. Record-writer stacks account for 32.9% and 34.0%; GPU waits account for 3.2% in the Metal trace. These are inclusive sample counts: nested categories overlap, and they are neither wall-clock stage durations nor GPU utilization. The source snapshot owner already uses APFS file cloning when available; this change does not remove source reconstruction or verification.

## Complete-cohort results

The [protocol](PROTOCOL.md) preceded candidate execution. Every run retains all **24,673 cells, 15,706 features and 14,184,532 records**. Each of twelve final old/new comparisons matches all four complete bundle files against the earlier qualified reference, and both final backends pass native reconstruction. The initial candidate's twelve comparisons also passed and remain archived. A subsequent safety review added a capacity guard after failed flushes; the table below is the final guarded owner, not the initial candidate.

| Owner / backend | Three elapsed times (s) | Median (s) | Mean (s) | Peak RSS range (MiB) |
| --- | --- | ---: | ---: | ---: |
| Prior CPU | 1.065, 1.070, 1.070 | 1.070 | 1.068 | 75.50–75.64 |
| Final CPU | 1.397, 0.906, 0.918 | 0.918 | 1.074 | 74.86–75.94 |
| Prior Metal | 1.070, 1.067, 1.062 | 1.067 | 1.066 | 86.83–92.05 |
| Final Metal | 0.912, 0.912, 0.913 | 0.912 | 0.912 | 86.28–86.38 |

Median elapsed time decreases by 14.2% for CPU and 14.5% for Metal. However, the first final CPU call is slower, and CPU's three-run mean is slightly worse. The initial candidate showed the same pattern: CPU 1.640/0.906/0.915 s and Metal 0.934/0.913/0.907 s. No run was excluded or retimed. This is one cohort on a shared physical M4 Pro desktop, using fixed rotated order and per-process source verification, publication and hashing. System caches were already exercised; cold-cache performance and the cause of the first-call penalty are not established. These observations do not establish general CPU acceleration or a meaningful Metal-over-CPU advantage.

## Safety and product checks

The actual writer passes release and AddressSanitizer checks for zero records, one record, 65,535/65,536/65,537 records and 131,075 records. An independent Python `<IIQ` encoder matches every byte from both builds. File and hash-only modes agree. Four out-of-range coordinates and overwrite reject in every case. A closed-output flush fails and the next append rejects before it can exceed buffer capacity. AddressSanitizer reported no error in these bounded Swift-owner checks; this is not full-application memory-safety qualification or sanitized C++ qualification.

All thirteen actual product single-cell CLI checks pass under the final executable identity: complete count import, default/explicit CPU identity, Metal output, three replays and malformed/overwrite/old-implementation rejection. The scoped build uses real owners and router via `Tools/Omics/H5AD/build.sh <output> --with-cli`; no generated substitute implements product behavior. It does not qualify the full application build.

CPU values retain SHA-256 `73d208729a6f63e355c61eedd2623eaf4b1547fd6196d6b23476976bda7585ce`; Metal values retain `3ebc3de24f974761afbc132c055cecbb9bc0d026f53c3563239fbaef009b3ee8`. Thus the earlier complete [Scanpy/NumPy precision checks](../README.md) remain exact numerical references for these payloads. Historical zero-tag fixture receipts are externally bound to actual executable/source hashes; the product CLI separately binds its real executable and OS fingerprint.

## Reproduction and evidence

The [archive manifest](evidence/2026-09-11/manifest.json) checks stored and decoded bytes. Its `source-map.json` maps every captured source file to deduplicated archive content. It retains raw profiles, both candidate comparisons, source/binary identities, release/ASAN and CLI checks, build warnings, receipts and cleanup mappings. Executables and large cohort payloads are retained externally by exact identity; source manifests for both candidates were reconciled, with the initial writer source archived separately.

- `Main.swift` / `run.py`: sample the retained owner and verify all sixteen receipts.
- `summarize.py`: extract the declared inclusive sample categories.
- `compare.py`: run the fixed full-cohort comparisons and native replays; `NUMIVIVO_PROFILE_ROOT` selects the final run directory.
- `WriterCheck.swift` / `check_writer.py`: exercise actual compiled owners and independently verify record encoding.
- `check_cli.py`: run the actual product router; `report.py` retains every timing and summarizes both candidates.

Generated duplicate values and CLI source copies were removed only after exact hashes matched retained files and `lsof` found no open handles. Receipts, logs and exact restoration paths remain. Restore the mapped payloads before replaying those bundles. Earlier Norman relocation and its verified compressed backup remain as documented in the [parent report](../README.md).

This is storage and numerical qualification. It adds no held-out biological outcome, does not qualify million-cell GPU execution, and does not close PCA/neighbors/model-fitting acceleration. [Biological prediction evidence](../../../../../Documentation/BiologicalPrediction.md) and the [full roadmap](../../../../../Documentation/SingleCellInteroperability.md) remain separate acceptance criteria.
