# Wider HNSW query: higher anchor recall, biological failure remains

Increasing query width from 128 to **512** in the fixed full-Gaussian MNN
candidate does not resolve Kang preservation. Megakaryocyte recall is **0.644707**
against original PCA's **0.750000**, a 0.105293 loss exceeding the unchanged 0.05
allowance. Condition balanced accuracy is **0.953524** against **0.973819**,
exceeding the unchanged 0.02 loss allowance. No candidate is promoted.

| Complete original cohort | Cells | Exact-anchor precision / recall | Fit seconds | Measured biological gates |
| --- | ---: | ---: | ---: | --- |
| Hagai | 13,863 | 100.000% / 99.996% | 3.816 | Pass |
| Kang | 24,673 | 99.999% / 99.984% | 15.723 | Condition and per-type recall fail |
| Ding | 44,031 | 99.971% / 99.922% | 32.988 | Pass |

All **82,567 cells** participate. Four original Kang treatment-classifier strata
remain unavailable, not passing; Ding retains partial source labels. The prior
width-128 candidate already failed Kang megakaryocyte recall. Its failure is
retained alongside this result. Closer approximate matching does not imply
monotonic improvement in preservation through panorama assembly and smoothing.
This test does not establish the causal contribution of any particular anchor.

## Fixed intervention and checks

The predeclared intervention changes only query ef in the native eligible-prefix/
suffix HNSW bridge. M=16, construction ef=200, seed=7, normalized Manhattan
matching, k=20, panorama assembly, sigma=15 and the qualified full-anchor Gaussian
library remain unchanged. Source cohort identities, original baseline values,
evaluation folds, programs and thresholds are retained. All cohorts were already
inspected development data; there is no independent validation claim.

Actual native fitting and replay complete on the physical Mac mini. Every replay
array is exact. The independent checker recalculates every selected matching
distance and full Gaussian update, reconstructs final coordinates and compares
all anchors with the exact native reference. Numerical checks pass on all three
cohorts. The biological evaluator runs in a separate pinned Python environment;
versions and executable/header hashes are retained. Fit times exclude witness
serialization/hashing and are not end-to-end application or scverse speedups.
No million-cell MNN run or production Swift integration is qualified here.

The post-fit numerical checker reports whether matching changed from the prior
candidate instead of asserting it unchanged. That assertion belonged to the
previous smoothing-only experiment. The exact fit-time source snapshot is
retained separately under `fit-source`, including the original checker, so this
diagnostic adjustment does not relabel the frozen fitting source.

## Retained evidence and next decision

`python3 verify.py` checks every archive member, fit-source/library identities,
output/array bindings and complete cohort gate accounting. The archive includes
native bridge/source headers, library, protocol, commands, environment versions,
all JSON results, evaluators and logs. [manifest.json](manifest.json) binds the
complete NPZ witnesses retained on the Mac mini, including all per-step Gaussian
arrays and biological neighbor outputs. Original source cohort artifacts remain
separate dependencies; the archive does not contain those raw datasets or all
large arrays. [summary.json](summary.json) retains every gate.

Restore those dependencies and external arrays to run the archived checkers.
For fresh fitting, use archived `FullGaussianMNN/run.py` with the recorded source,
manifest, native and Gaussian libraries, `PROTOCOL_WIDE.md`, and a new output
path. Then use the post-fit `check_numerics.py` and `LocalMNN/evaluate.py` with the
recorded evaluation spec; adapt absolute paths without changing input hashes.

Further width tuning on Kang would be development reuse. The remaining problem
requires investigating how small anchor changes affect assembly and cell-type
boundaries, while retaining all original biological gates. Neither near-exact
anchor recall nor the faster approximate matcher is sufficient for promotion.

The subsequent [exact-anchor assembly control](../ExactAnchorControl/README.md) reproduces native coordinates within 1.3e-13 on all three cohorts. Kang's 26 missing anchors involve 22 distinct reference-labeled megakaryocytes; aggregate recall obscures this concentration.
