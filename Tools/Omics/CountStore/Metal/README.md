# Explicit Metal sparse-count normalization

The native count-store owner now offers `metal-fp32` alongside the unchanged default `cpu-fp64`. Complete original Kang qualification retains **24,673 cells, 15,706 features and all 14,184,532 source-major records**. This implements an optional Apple GPU transform; it does not qualify downstream biological prediction or promote a production backend.

```sh
numivivo singlecell-count-store-normalize <store> --target 10000 --output <new-directory> --backend metal-fp32
numivivo singlecell-count-store-normalize-verify <normalized> --store <store>
```

Omitting `--backend` retains the CPU path. The actual product single-cell router was compiled and exercised with its real executable-fingerprint owner through `Tools/Omics/H5AD/build.sh <output> --with-cli`. This scoped executable is named `numivivo-omics`; it is not a full NumiVivo application build. The H5AD harness separately checks original zero-tagged historical fixtures, with actual binary and compiled-source hashes recorded in evidence.

A final default-mode wrapper check exposed Bash 3 empty-array handling and was repaired. The rebuilt executable has a different hash despite identical recorded source files; both complete CPU and Metal bundles replay exactly under that build. The failed binary-identity expectation, successful replay and both identities are retained. Timings below remain bound to the original benchmark executable.

## Precision, storage and verification

Raw counts remain exact UInt64. The Metal path derives row scales in FP64, converts scales and counts to FP32, applies a compensated FP32 log1p with fast math disabled, and exactly widens results into the existing FP64 record container. Its optional execution receipt records the precision profile, kernel hash, physical device/registry ID and 262,144-entry batch size. CPU receipts omit that record and remain byte-exact to the retained pre-change owner.

Input/output GPU buffers are reused after each completed command: 4 MiB input and 1 MiB output, plus one FP32 scale per cell. Count reads retain the existing 16 MiB mapping window and output uses the existing 1 MiB hashed writer. There is no dense cells-by-genes allocation. Metadata and quality arrays remain resident. Targets are restricted to 1...1e9 for Metal; the CPU target range is unchanged. Missing/nonphysical GPU support fails explicitly. No CPU fallback or receipt is published after a failed GPU command.

Native verification reconstructs the source and reruns the recorded backend; device/profile differences reject verification. This intentionally pins GPU reconstruction to the recorded device. It does not assert portability of exact results across arbitrary GPUs or operating systems.

## Complete-cohort numerical results

| Comparison | Records | Maximum absolute error | Violations |
| --- | ---: | ---: | ---: |
| Mac mini CPU vs FP64 Scanpy | 14,184,532 | 1.78e-15 | 0 |
| Mac mini Metal vs FP64 Scanpy | 14,184,532 | 1.24e-06 | 0 |
| Local M4 Metal vs FP64 NumPy | 14,184,532 | 1.24e-06 | 0 |

Coordinates and metadata match exactly. Metal uses the predeclared tolerance `3e-6 + 3e-6 * abs(reference)`; CPU uses `1e-14` absolute. All Metal values are positive, finite and exact FP32-to-FP64 widenings. All three same-device publications match byte for byte; both backends pass native reconstruction. The old default owner and new CPU owner produce identical complete bundles.

Controlled physical-GPU checks cover six full/partial-batch runs, UInt64-scale totals and three normalization targets. Nine rejection checks cover invalid targets, overwrite and altered device/profile/backend records. Both inspected hosts have physical Apple GPUs, so the unavailable-GPU branch was not exercised. The primary numerical evidence is the complete experimental cohort, not these controlled cases.

## End-to-end owner timings

These physical M4 Pro desktop measurements include source snapshot verification and reconstruction, metadata, per-process kernel setup/conversion, normalized output and hashing. CPU and Metal were alternated. Earlier controlled checks warmed system caches; no cold-cache or isolated-kernel claim is made.

| Run | Elapsed seconds | Peak RSS MiB |
| --- | ---: | ---: |
| cpu-before | 1.500 | 76.53 |
| cpu-0 | 1.394 | 81.14 |
| metal-0 | 1.074 | 87.02 |
| cpu-1 | 1.039 | 75.58 |
| metal-1 | 1.050 | 86.86 |
| cpu-2 | 1.043 | 75.56 |
| metal-2 | 1.048 | 92.03 |
| cpu-0-verify | 1.115 | 75.77 |
| metal-0-verify | 1.134 | 87.05 |

Three-run medians are **1.043 s CPU** and **1.050 s Metal** (Metal/CPU elapsed ratio 1.007). Metal is slower in this complete-pipeline observation. The sample is one dataset on a shared desktop; no general performance advantage, default change or biological acceptance follows.

The independent Scanpy transform took 1.126 s after 0.089 s of H5AD loading on another host. It omits the native publication/reconstruction boundary and is therefore a descriptive reference timing, not a matched speed ratio.

## Evidence and restoration

[The frozen protocol](PROTOCOL.md) and [archive manifest](evidence/2026-09-11/manifest.json) retain source/binary identities, every run, independent checks and failure coverage. `run_complete.py`, `check_complete.py`, `check_rejections.py` and the two Swift drivers reproduce the checks using the retained source paths and owners.

Before qualification, two identical historical Norman normalization payloads (5,785,321,936 bytes each) were preserved in one local compressed copy, fully decoded and checked against both receipts and files. Only then were the two remote payloads removed after empty open-handle checks, recovering 5,785,112,576 measured available bytes. The archive records its exact backup path, compressed and decoded hashes; old receipts and metadata remain. Restore either payload from the verified gzip copy before replaying its historical normalized bundle. This storage relocation changes no old numerical evidence.

PCA, neighbors, model fitting, million-cell GPU execution, matched end-to-end scverse comparisons and downstream biological preservation remain separate work. See the [complete roadmap](../../../../Documentation/SingleCellInteroperability.md).
