# Integration scratch lifetime — 10 September 2026

Integration publication now removes original-score, normalized-score and distance
scratch immediately after solving. It then serializes and fingerprints each
output matrix before removing that matrix's scratch. Unconditional, idempotent
owner cleanup still handles failures and cancellation. Solver arithmetic,
artifact bytes, report bounds and transaction boundaries are unchanged.

## Qualification

The complete native release build and 16 tests in four suites passed, including
successful scratch consumption with exact output bits, write failure preserving
scratch for the owner, cancellation cleanup, and prior batch/trajectory checks.
All 533 production source hashes bind executable SHA256:
`e360645362ac0707cb497eaa64277e5bf792b28d79bddea707421fdcd043495d`.

Fresh complete Kang/Hagai runs (24,673/13,863 cells; seeds 7, 19 and 41) reproduce
all previously qualified matrix payloads and numerical diagnostics. All six
native replays pass. Another 60 ridge lifecycle commands pass with 16 expected
rejections; 45 MNN commands pass with 11 expected rejections and the retained
independent Scanorama numerical checks. Mac mini execution used macOS build
25G72; fresh laptop fixtures retain their actual build 25G5028f identity.

## Full-shape lifetime measurement

A synthetic storage driver uses the unchanged production matrix APIs with
1,612,594 rows, four 20-column matrices and two 100-column matrices. Both
schedules write all three complete outputs; each output is independently read
back and hashed. Every output hash agrees between schedules.

| Quantity | Old schedule | Early retirement |
|---|---:|---:|
| Calculated peak scratch plus output extents | 8,566,099,328 B | 5,160,300,800 B |
| Sampled peak logical file extents | 8,566,099,328 B | 5,160,300,800 B |
| Sampled peak allocated blocks × 512 | 8,235,393,024 B | 4,748,922,880 B |

The logical peak falls by 3,405,798,528 bytes. The 25 ms sampler observes lower
bounds on transient peaks; matching the layout calculation provides a separate
check for logical extents. Allocated-block accounting is different: APFS can
report zero allocated blocks while dirty mapped scratch still occupies memory.
These measurements exclude parent/PCA reconstruction, metadata, unrelated
filesystem usage and swap. They are not whole-run disk or RSS bounds.

The final driver completed successfully in 24.79 seconds and removed all its
scratch and synthetic outputs. This is a storage lifecycle measurement, not a
full integration solver, biological or controlled performance qualification.

The first observer failed when native cleanup removed a directory during
iteration. Its native result reported matching outputs, but the wrapper did not
observe native exit status. That incomplete observation, its sources and logs
are retained separately. The corrected observer handles directory disappearance;
its second run supplies the successful terminal status and sampled peaks above.

## Next full-cohort gate

Full HIRISA integration still requires fresh PCA ancestors under this runtime,
complete native publication/replay, independent numerical reconstruction and
biological-preservation evaluation. The known Kang NK-cell preservation failure
is unchanged. The predeclared operational baseline retains all original cells,
uses 100 clusters, donor covariates, seed 7 and at most 10 iterations. It is not
an assertion of general integration competitiveness or prospective prediction.

## Archive

The [manifest](evidence/2026-09-10-integration-scratch/manifest.json) retains native
checks, all small-cohort payloads, lifecycle fixtures, both lifetime observation
attempts and runtime/source identities. Existing compressed matrix payloads are
reused. Large original H5AD inputs and executables remain external with hashes.

The manifest contains 1,510 members and 217,656,187 logical compressed bytes.
Of these, 217,407,548 bytes reused verified existing payloads when collected;
only 329 new compressed members needed transfer. Manifest SHA256:
`365fa7694b509b66ca90aea61250edc9eb5b43b94b5d8768937a3e06e1f3e0c2`.

Verify with `python verify_archive.py evidence/2026-09-10-integration-scratch`.
