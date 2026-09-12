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
