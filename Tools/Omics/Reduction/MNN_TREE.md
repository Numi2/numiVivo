# Exact-tree MNN trial: identical results, slower execution

An exact spatial-tree matcher preserves every original MNN score and anchor byte
on the complete Hagai, Kang and Ding cohorts, but is slower in all three measured
runs. **The prototype was not promoted.** The production API, matcher, numerical
parameters and defaults remain unchanged. Its full source overlay and results are
retained so this approach need not be repeated without a new reason.

The remaining production change repairs the scoped H5AD build's explicit source
list: `VivoBufferedCountRecords` was missing although the existing clustering
owner depends on it. The restored production build passes on the physical Mac
mini. No numerical implementation changes accompany that build repair.

## What was tested

The candidate replaced exhaustive cross-level matching with one balanced spatial
tree over the same normalized PCA coordinates. Queries search strictly lower or
higher donor/batch levels, retaining the existing Manhattan metric and grouped
source-row tie order. Bounding boxes prune only when their lower distance bound
is strictly greater than the worst retained distance; equal-distance candidates
remain eligible. Every cell participates, including duplicate/zero directions.
The Gaussian correction still uses all original anchors and the same median-norm
scale, alignment selection, panorama assembly and reduction order.

The native bridge tracks construction, bounding-box and point-distance coordinate
terms against a common work limit. Its allocated index estimate and the bridged
level array are included in latent admission. Exact search remains worst-case
quadratic; this trial does not establish million-cell MNN integration.

All inputs were checked against the [published native MNN cohort manifest](evidence/2026-09-10-native-mnn/manifest.json)
before execution and rehashed afterward. One scan and two tree executions used
the same statically linked production numerical owner on a physical M4 Pro.
These timings exclude input decoding/output publication from the owner interval;
full command logs and OS peak RSS are retained separately. They are not repeated
end-to-end CLI timing distributions or GPU measurements. The separate original
HIRISA clustering baseline remained active; no competing job was cancelled.

| Complete cohort | Cells | Scan owner seconds | Tree owner seconds | Tree replay seconds |
| --- | ---: | ---: | ---: | ---: |
| Hagai | 13,863 | 4.509 | 5.333 | 5.351 |
| Kang | 24,673 | 31.124 | 34.979 | 35.042 |
| Ding | 44,031 | 68.607 | 92.421 | 92.920 |

The tree adds work for bounding boxes and performs separate directional queries,
whereas exhaustive matching reuses each pair distance in both directions.
Point-distance work decreases for Kang but increases for Hagai and Ding.
The Gaussian kernel's scalar-term counts are unchanged: 4.934, 41.258 and 88.536
billion respectively. On these cohorts, replacing the matcher alone is therefore
not sufficient to establish faster or scalable integration.

## Numerical evidence and retained limits

All nine executions match the original complete corrected-coordinate and anchor
files **byte for byte**, including every anchor distance. All alignments, assembly
orders, per-step corrections, zero-weight counts and panoramas match exactly.
The current exhaustive mode also reproduces the original report exactly. Tree
report replay is exact, including work counters.

Twenty-four independent NumPy exhaustive-query cases cover 1, 2, 20 and 64
dimensions, ties, mixed source-row/level order, and 1, 20 and 100 neighbors. Every
lower/upper neighbor is checked. Work-boundary, cancellation, insufficient-output
capacity and nonfinite-input checks pass. Four native Swift tests pass, including
the original independent Scanorama fixture, uniform-scale behavior, admission
failures, legacy serialization and tree-versus-scan coordinate bits. Address and
undefined-behavior sanitizers pass local balanced-leaf boundary and error-path
checks; these are numerical/memory checks, not biological validation.

Identical score bytes preserve the existing [cohort metrics and limitations](MNN_INTEGRATION.md).
Kang's unavailable classifier strata, Ding's partial labels, inspected-study
selection and the absence of independent biological validation remain unchanged.
HIRISA's 1.61-million-cell integration and annotation failures are a separate
experiment; this tree trial does not resolve them.

The first harness execution exceeded the PCA reader's tile limit before running
correction. It was changed to 2,048-row reads. Initial scoped-build and standalone
Testing framework/plugin/fixture-path failures are recorded with their fixes.
No cohort, numerical parameter or acceptance margin changed between attempts.

## Reproduction and next algorithm work

The [archive](evidence/2026-09-11-mnn-tree) contains 75 deterministic gzip members,
116,871 stored bytes. Its manifest SHA-256 is
`5deef18182f4674f446659536669442e7e16f4a4541e57b6fdb927e04404fbf4`.
Verify stored and decoded hashes with:

```sh
python Tools/Omics/Benchmarks/HIRISA/verify_archive.py Tools/Omics/Reduction/evidence/2026-09-11-mnn-tree
```

Restore the archived `prototype/` overlay into a separate clean checkout of base
`435531e369f72f343663bfeeb442f2f71340ceb7`; verify source-manifest.json first.
Build with its H5AD build script. Link MNNTreeMain.swift against that scoped
library plus OmicsHNSW.o and OmicsMNNTree.o; the production source bindings and
compiler are recorded in environment.json. The tested harness SHA-256 is
`e96d590308192551ca2de03d5e7fd047ffd464f95a15a6940896c9ab92a56619`.
Run the archived run_mnn_tree.py with the original cohort root, prior-artifacts.json
and a new output directory. The full coordinate and anchor outputs are omitted
from this small archive because each exactly matches a hash-bound original MNN
payload; the manifest gives every restoration mapping. The prototype is research
evidence and is intentionally absent from the production API.

The next scalability design must address both neighbor search and Gaussian
correction cost. Any approximate matching or restricted-kernel method needs an
explicitly frozen numerical contract, all original cells, anchor/coordinate error
measurements, and the unchanged biological-preservation diagnostics before
promotion. This result does not justify an approximate option without those checks.
