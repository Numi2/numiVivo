# Where rigid integration loses annotation recall

Recounting the existing full-HIRISA evaluation places **12 of 31 failures in
B intermediate (7/8 eligible) and CD14 Mono (5/8 eligible)** comparisons.
CD8 TCM fails 4/4, B memory 3/8, and eight other labels account for the remaining
12 failures. These are label-by-stratum comparisons, not independent experiments
or counts of misclassified cells. Rare status belongs to each comparison; a
label can have both rare and non-rare comparisons.

[diagnosis.json](diagnosis.json) retains every label, eligible and unsupported
counts, insensitive supported counts, failed stratum IDs and maximum recall
losses. All 471 comparisons remain represented: 146 are eligible, 31 fail,
and 9 of the failures are rare. Unsupported comparisons are not passes.

Run `python3 diagnose.py > regenerated.json` and compare the parsed JSON with
`diagnosis.json`. The script verifies the parent archive hash and the exact
annotation evaluation member before aggregating. It neither reruns the model
nor independently requalifies the original scorer. The parent experiment's
numerical checks and biological failures remain unchanged.

## Next integration experiment

Inspect the retained donor-fold confusion for B intermediate and CD14 Mono
alongside the rare failures to distinguish displacement across related labels
from broad loss of separability. This is a diagnostic use of author labels,
not permission to fit a correction to the held-out labels. Keep all original
annotation and program comparisons in any subsequent candidate evaluation.

Before fitting another candidate, freeze its training-only selection rule,
source coordinates, control groups, donor exclusions and complete acceptance
rules. Compare with original PCA and both existing corrections on the same
folds; include the mixed-program failures as well as annotation recall. A
candidate must address cross-donor geometry: preserving distances separately
inside each donor has already proved insufficient. Do not promote a method
based on fewer failures, improved donor mixing or a passing subset. These
inspected HIRISA data remain development evidence; any claimed general
biological preservation requires a separately frozen validation cohort.

## Donor-fold destination trace

[confusion.json](confusion.json) preserves all 471 comparisons and their fold
rows, with baseline and candidate destination counts. `python3 confusion.py`
regenerates it from the hash-checked parent archive. It checks all 114 matched
folds, integer nonnegative confusion counts, source query/training denominators,
recalls and count conservation. No annotation model is fitted.

Across sufficient folds of all eight eligible B intermediate comparisons,
correct assignments fall from 46,829 to 39,668. Assignments to B naive increase
by 4,415 and to B memory by 2,758. Recall falls in 33/39 sufficient folds,
with losses present in every donor. Across the eight eligible CD14 Mono
comparisons, correct assignments fall from 238,135 to 236,722; assignments to
CD16 Mono increase by 4,568, partly offset by reduced assignment to other labels.
Recall falls in 19/39 sufficient folds, also spanning every donor.

These pooled counts describe query appearances across evaluation strata, not
independent biological replications. They do not replace the original equal-fold
recall gates. Confusion matrices provide marginal counts: they cannot identify
which individual cells changed prediction between the two methods, and the
changes do not establish a causal mechanism. Author labels remain references
for this diagnostic, not authoritative automatic assignments.

The observed loss concentrates in distinctions within the B-cell and monocyte
labels, motivating a check of local cross-donor neighborhood alignment before
another global correction. A label-guided fix or selective relabeling would not
establish preservation. Retain the complete original gates and require a fresh
cohort before claiming generalization.
