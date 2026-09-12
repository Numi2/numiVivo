# Exact-anchor control isolates the matching divergence

The same full-Gaussian research assembly reproduces native exact-MNN coordinates
when supplied the native exact anchor sets. All **82,567 original cells** are
retained, and assembly order matches on all three cohorts.

| Cohort | Cells | Exact anchors | Maximum coordinate error |
| --- | ---: | ---: | ---: |
| Hagai | 13,863 | 49,893 | 1.51e-14 |
| Kang | 24,673 | 159,513 | 1.30e-13 |
| Ding | 44,031 | 131,282 | 2.33e-14 |

This executes every full-Gaussian correction step using original PCA and the
qualified native kernel. It does not rerun exact anchor discovery: the frozen
native anchors are intentionally the controlled input. No cell labels enter
assembly. This isolates the larger approximate-coordinate divergence to matching
for these fixed cohorts, within the measured numerical tolerance. It is not a
new scalable matcher, biological qualification or product promotion.

## Changed anchors and rare-cell concentration

In the preceding width-512 candidate, Kang misses **26** exact anchors and adds
two. The missing anchors involve **28 distinct cells**, of which **22** carry the
reference megakaryocyte label, five CD14+ Monocytes and one FCGR3A+ Monocyte.
There are 52 missing-anchor endpoint occurrences, including 32 megakaryocyte
occurrences; repeated endpoints are not independent cells. The added anchors
involve two megakaryocytes and one CD4 T cell.

All changed anchors from all cohorts remain counted: Hagai misses two and adds
none; Ding misses 103 and adds 38. Unreported labels remain explicit, including
all Hagai endpoints and unlabeled Ding cells. These are retrospective source-label
diagnostics, not authoritative annotations or a causal intervention on individual
anchors. The full exact-anchor control does not identify which missing or added
anchor causes a particular classification change.

Aggregate anchor recall hides a concentrated rare-cell matching problem. The
next matching method should evaluate coverage and stability in sparse local
regions, using geometry rather than held-out labels to decide computational
effort. Any adaptive rule must be frozen before fitting and retain all original
biological gates and missing strata. Selecting corrections using the inspected
Kang labels would require new independent validation and would not establish
unsupervised preservation.

## Evidence and reproduction

[results.json](results.json) retains every full-coordinate comparison;
[changed-anchors.json](changed-anchors.json) retains all population counts and
source bindings. `python3 verify.py` verifies the archive, frozen source/protocol,
complete cohort accounting and endpoint count conservation. It does not rerun
numerical assembly or discover anchors.

The archive retains the predeclared protocol, executed fitting/driver scripts,
per-cohort reports, input/library hashes, post-run endpoint diagnostic and logs.
The complete per-step NPZ witnesses are retained on the Mac mini and bound by
[manifest.json](manifest.json). The original source cohorts, wider-search candidate
and qualified Gaussian library remain external dependencies. Restore those exact
inputs and use the archived `run.py` in a fresh output directory, adjusting
absolute paths without changing identities. `changed_anchors.py` separately
recounts the preceding approximate candidate's changed endpoints.
