# Protected-stratum regression: measured tradeoffs, no promotion

One fixed-membership development experiment on every original **1,612,594 HIRISA
cell** does not repair the published annotation-retention failure. Adding separate
preparation/treatment intercepts leaves **37/146** supported, control-sensitive
label/stratum comparisons outside the unchanged recall-loss margins. Four failures
resolve and four new failures appear. No native integration default or solver has
been changed on this evidence.

The [protocol](PROTECTED_REGRESSION_PROTOCOL.md) was declared after inspecting the
[original annotation result](ANNOTATION_RETENTION.md). It preserves all twenty PCs,
the original seed-7 one hundred soft-cluster memberships, five donors, ridge one,
and the original active-donor criterion. The single alternative estimates donor
effects conditional on 23 original experimental preparation/treatment strata;
only donor effects are subtracted. Author labels, confidence and RNA program or
response outcomes are not used to fit the correction. Frozen memberships isolate
the regression change, but this is **not an iterative native integration refit**.

## Results on unchanged diagnostics

| Endpoint | Original native | Protected regression |
| --- | ---: | ---: |
| Annotation failures / supported, sensitive comparisons | 37/146 | 37/146 |
| Rare annotation failures / supported, sensitive rare comparisons | 10/29 | 8/29 |
| Mean recall over the 146 eligible comparisons | 0.650832 | 0.652043 |
| Coarse response contrasts meeting all declared margins | 16/16 | 16/16 |
| Within-library program preservation among sensitive comparisons | 28/28 | 28/28 |
| Within-library program controls lacking sensitivity | 4/32 | 4/32 |
| Original mixed-objective program failures among sensitive comparisons | 3/14 | 3/14 |
| Mean conditional donor scatter fraction over 23 groups | 0.020288 | 0.016428 |
| Mean relative PC response drift over 16 contrasts | 0.038105 | 0.036461 |

Means weight comparisons or groups equally, not individual cells. Original PCA's
mean conditional donor scatter fraction is 0.079429. Reduced donor-associated
variance cannot identify technical batch removal separately from donor biology.
The complete annotation and program gates remain false. All 471 annotation
comparisons, including 325 without sufficient support or control sensitivity,
remain in the evidence; rare means below one percent within the original stratum.

Among the 146 eligible annotation comparisons, recall improves in 73, worsens in
58 and is unchanged in 15. The average improvement is **0.121 percentage points**;
individual changes range from **−5.889 to +4.557 percentage points**. Resolved
failures are HSPC in Monocyte IFNb and IFNg, CD8 TEM in Tcell IFNa, and NK in Tcell
none. New failures are CD4 TCM in NK none and CD16 Mono, CD4 Naive and CD8 TEM in
PBMC Fresh. The lower rare-failure count therefore does not imply a general repair.

This was one hypothesis on inspected data, with no strength sweep, altered margins
or relabeling. It supplies evidence about one component of correction and its
tradeoffs. It neither proves that protected covariates are generally ineffective
nor supports production promotion. Subsequent work needs a separately specified
algorithm and the full preservation evaluation, including neighborhood structure
and independent contexts. It must retain these failures and insufficient strata.

## Numerical checks and reproduction

The original native coordinates were reconstructed across every cell with maximum
absolute discrepancy **1.8833e-11**. Every protected cluster regression was checked
against independent augmented least squares on weighted stratum/donor means;
the resulting full coordinates differed by at most **7.2733e-12**. This grouped
regression oracle is separate from the annotation oracle: baseline, reconstructed
native and protected annotation fits were each checked by per-cell weighted SVD
on all 114 donor-held-out folds, with exact query confusion. Reconstructed native
confusion and metrics exactly match the earlier publication. All 79 matched
within-library program folds for both programs pass the independent per-cell
oracle. Three focused regression tests cover the original single-intercept case,
complete confounding and invariance of donor effects to stratum offsets.

Use the qualified NumPy/h5py environment and the study-relative inputs bound by
the execution freeze. The output directory must be new, with protocol.md copied
from the declared protocol. Run in order:

```sh
python fit_protected_regression.py --study STUDY --out STUDY/protected-regression
python evaluate_protected_regression.py --study STUDY --root STUDY/protected-regression --python PYTHON
python summarize_protected_regression.py --study STUDY --root STUDY/protected-regression
python archive_protected_regression.py --study STUDY --root STUDY/protected-regression --out ARCHIVE
python verify_archive.py ARCHIVE
```

The [evidence archive](evidence/2026-09-11-protected-regression) preserves models,
freezes, all results, original-row metadata, evaluator sources and independent
checks: 74 members, 7,309,063 stored bytes, manifest SHA-256
`b05bf49d703b114687c44419e6cb9c36f35e9fc88d2152d3e222fe059c63e40f`.
Both full coordinate payloads remain external, each 516,030,080 bytes:

- Reconstructed native SHA-256: `b5ac524ab722aa442431d95b1958bb11423dedf2003c712b8782892916692e5c`.
- Protected SHA-256: `35ab3be48ebdd0dbcf9ecea9aee6f2097a528dc35af62ce8741604d00a6f07b3`.
- Models SHA-256: `625867fa7b62e5cac794e0b393a52b74e7ac32653c8675f4a243379cd7b2b115`.

The archive records original PCA/membership/source bindings and external payload
restoration paths. Initial summary/archive invocations failed before archive
creation because system Python lacked h5py and the dependent summary was absent;
a subsequent archive attempt also caught frozen scripts not yet copied into the
isolated repository. Using the qualified environment and copying the identical
scripts resolved these failures without rerunning or changing scientific outputs.
All attempts are retained in publication-attempts.json.
This transductive author-label diagnostic is not prospective biological prediction,
independent validation, clinical efficacy or a variant-to-phenotype qualification.
