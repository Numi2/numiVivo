# First sparse-normalization Metal qualification

Declared on 2026-09-11 before running the candidate on experimental counts.
Base: `9cb829b8a800e3e89e0437d912d4e500234a2148`.

Run the complete original Kang analytical count store: 24,673 cells, 15,706
features and 14,184,532 source-major count records. Preserve the qualified RNA
mapping, source order, every cell/feature identity, raw count and metadata.
Target total is 10,000. No dense cell-by-gene array is permitted.

The existing scalar FP64 path remains the default. The explicit Metal path
computes row scales in FP64 then converts them and counts to FP32; compensated
FP32 log1p values are widened exactly into the existing FP64 container. Retain
the numerical profile, kernel hash, physical device and fixed batch size in its
receipt. Require a physical Apple GPU, no silent CPU fallback, target 1...1e9,
positive finite outputs and successful command completion before publication.
Historical CPU receipts omit the new optional execution record.

Compare every output record against independent NumPy/Scanpy normalization,
with exact coordinates and an absolute plus relative tolerance of
`3e-6 + 3e-6 * abs(reference)` for Metal. Require FP64 CPU agreement within
`1e-14` absolute error. Report maximum errors and every violation, and verify
source reconstruction and exact same-device native replay for both paths.
These numeric margins are not a biological preservation claim.

Time three complete CPU and three complete Metal publications in alternating
order, including source snapshot verification/reconstruction, metadata, kernel
creation, conversion, output writing and hashing. Retain the first Metal publication
separately. Every run creates a fresh process and kernel owner, but earlier
controlled kernel checks warm system caches; this is not a cold-cache benchmark. Record process elapsed time and peak RSS. Record an independent
full-data Scanpy normalization time and memory with its different I/O boundary
explicit; do not present it as a matched end-to-end comparison. Do not infer a
speedup from kernel-only timings. Preserve all results if Metal is slower.
No default or downstream inference behavior is promoted from this experiment.

Controlled checks cover a partial final GPU batch, low count/large total ratios,
large UInt64 counts, invalid targets, unavailable GPU, modified receipt/device
identity and overwrite protection. They support arithmetic and artifact safety;
the complete original cohort is the operational benchmark. Keep external active
GPU jobs untouched and check workload ownership immediately before timings.
