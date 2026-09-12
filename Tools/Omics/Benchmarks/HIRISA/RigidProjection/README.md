# Full-HIRISA rigid projection: preservation still fails

A donor-wise rigid projection on all **1,612,594 cells and 20 PCs** reduces the
annotation failure count from 37 to **31/146**, but does not repair preservation.
Ten failures resolve and four new failures appear. Rare failures change from
10 to **9/29**. Mean recall across eligible comparisons is 0.650730, compared
with 0.650832 for the original native correction; fewer threshold failures do
not imply higher average recall. No production method or default is changed.

This extends the earlier three-cohort rigid-projection development method to
HIRISA. For each donor, fit a centered orthogonal Procrustes transformation from
original PCA to the frozen seed-7 native integrated coordinates, then restore
the target donor mean. Rotations/reflections are allowed; scaling is not. Author
labels and program outcomes do not enter fitting. All within-donor distances
are preserved algebraically, but different donor transformations can still
change cross-donor relationships and annotation recoverability.

| Unchanged endpoint | Rigid projection |
| --- | ---: |
| Annotation comparisons retained / eligible | 471 / 146 |
| Eligible annotation failures | 31 |
| Eligible rare failures | 9 / 29 |
| Coarse response contrasts passing | 16 / 16 |
| Sensitive within-library program comparisons passing | 28 / 28 |
| Within-library comparisons lacking control sensitivity | 4 / 32 |
| Sensitive mixed-objective program failures | 2 / 14 |
| Mean conditional donor scatter fraction | 0.028972 |
| Mean relative response drift | 0.050339 |

Conditional donor scatter averages the same 23 groups equally. Original PCA
has 0.079429, native correction 0.020288 and the prior protected regression
0.016428. Rigid projection leaves more donor-associated scatter than either
correction. These descriptive quantities cannot separate technical batch
removal from donor biology. All insufficient comparisons remain in the archive;
the complete annotation and program gates remain false.

## Numerical scope and reproduction

Inputs are hash-bound to the existing full-cohort annotation freeze. NumPy SVD
and SciPy polar decomposition agree for all five transforms. All corrected rows
are checked against separately applied polar transforms; the maximum relative
centered-norm discrepancy is 4.26e-15. This is matrix-method consistency plus an
isometry check, not an independent scientific reference. All 114 annotation
folds are independently checked by the existing per-cell SVD oracle, with exact
query confusion. The baseline reproduces prior confusion and metrics. The
within-library program oracle also passes across all original matched folds.

The protocol and inputs were frozen before fitting, and coordinates before
scoring. This is one Python development experiment on inspected data, not a new
native execution or independent prospective validation. A preserved donor's
geometry does not establish cross-donor biological preservation.

Run `python verify.py` to verify the archive and recompute the reported failure
counts. The archive retains the protocol, models, evaluators, results and failed
attempts. `fit.py STUDY OUT` and `evaluate-resume.py --study STUDY --root OUT
--python PYTHON` reproduce execution with the retained study-relative inputs;
use a fresh output directory, copy the protocol and evaluator dependencies, and
adjust archived absolute paths for another host. The manifest binds the external
516,030,080-byte coordinate payload, retained on the Mac mini. It is not embedded
in the 3.7 MB result archive.

The first dependency attempt selected an unsupported source build; SciPy 1.16.3
binary installation resolved it. The first evaluator stopped before scoring on
an output-directory assumption; only that path handling was corrected before
resuming, without refitting. Both failures are retained.
