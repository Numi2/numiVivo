# One full-cohort local MNN development candidate

Declared 2026-09-11 after the protected-regression and exact-tree failures.
This is development on inspected cohorts, not independent biological validation.
Do not tune the following settings or remove cohorts/classes after seeing results.

Retain every original Hagai, Kang and Ding cell, the exact original twenty-PC
payloads, original donor/batch assignments, k=20 mutual neighbors, sigma=15,
alignment threshold 0.1, global median-row-norm scale, panorama ordering and merge
multiplicity. No author cell labels, RNA programs or held-out outcomes enter fit.
Condition/covariate eligibility remains that of the already qualified original
cohorts. This candidate does not claim a prospective unseen-donor transform.

Use the pinned vendored hnswlib 0.8.0 for serial native FP64 neighbor queries:
M=16, construction ef=200, query ef=128, seed=7 for each independently constructed
index. Matching uses one index over all row-normalized original PCs, with strict
lower/upper covariate-level result filters. Sort returned candidates by Manhattan
distance and grouped original-row rank; this does not guarantee global exact ties.
Retain every mutual pair found. Search is approximate and is reported as such.

For each correction step, snapshot every selected anchor's source coordinate and
bias in original median-norm units. Build a separate squared-Euclidean index over
all anchor records, preserving duplicates. For every target cell, use at most the
64 returned nearest anchor records, or all records when fewer than 64 exist. Sort
returned records by squared distance and original anchor-record index. Apply the
same nonnegative Gaussian weights exp(-0.5*sigma*squaredDistance), normalized over
those selected records. Zero total weight preserves the target cell and is
reported. This is a deliberately local correction, not the original all-anchor
Gaussian kernel with an asserted approximation guarantee. No source cohort is
subsampled, and no anchor record is removed from the search index.

Admission: 2 million indexed rows per index, at most 64 components, at most 128
requested neighbors, and 500 million metric evaluations per index lifetime.
Retain index storage, all construction/query evaluation counts, all selected
kernel IDs/distances and all corrected coordinates. These bounds are not a
million-cell, out-of-core or GPU qualification. Failure preserves evidence.

Freeze source/input identities before fitting, then output identities before
reading biological metrics. Refit with the same fixed seed and verify every
coordinate/anchor/selected-kernel payload exactly. Independently recompute every
selected kernel distance and weight/update from original step snapshots and
check finite results. Compare all mutual anchors with the exact original native
anchor set; require global precision and recall >=0.95, and report every pair,
including changed alignments. Report full coordinate errors against the exact
original method without calling the local kernel an exact approximation.

Apply the unchanged complete-cohort 30-neighbor mixing, held-donor/experiment
classifiers, every-type recall, measured-program and within-stratum-program
margins. Preserve all insufficient/unavailable strata and partial source labels.
Measure timing on the same physical Mac mini; do not infer end-to-end release
speedup from a development driver. Even favorable numerical/biological results
only justify a separately implemented and qualified native method plus a full
HIRISA trial; they do not establish production or biological completion.
