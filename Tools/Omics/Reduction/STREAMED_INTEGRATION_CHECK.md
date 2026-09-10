# Streamed integration reconstruction

`check_streamed_integration.py` independently checks native categorical ridge
integration without allocating a resident cells-by-clusters or cells-by-genes
matrix. It reads at most 8,192 rows from each record file per batch. Per-cell
level indices, the bounded report, and small cluster/batch systems remain resident.

```sh
python check_streamed_integration.py --bundle INTEGRATION_DIRECTORY --out checks.json
```

Dependencies qualified here: NumPy 2.5.3 and ijson 3.4.0.post0. For a study with
predeclared metadata, add `--expected-metadata-sha SHA256`. The checker validates
all matrix coordinates and finite values, membership bounds and row sums,
per-cell covariates against original metadata, the final diversity objective,
and every corrected coordinate. Weighted sufficient statistics are accumulated
in row batches. Ridge coefficients come from an independent dense normal-equation
solve, rather than the native Schur-complement implementation. Every input is
hashed before/after numerical checking to reject changes during the check.

The checker covers fixed and expected-mass ridge, donor and batch covariates,
and fitted/query PCA inputs. MNN has its separate checker. This check establishes
final objective/correction consistency; it does not reproduce the optimization
trajectory, prove biological preservation or establish prospective prediction.
Native artifact verification/replay remains a separate required gate.

Membership row sums allow absolute error 1e-10. Objective and corrected-coordinate
comparisons use absolute 1e-8 plus relative 1e-9 to accommodate independent
summation and linear-solve order. These tolerances were frozen before the full
HIRISA integration output. They are numerical checks, not biological thresholds.

The [qualification archive](evidence/2026-09-10-streamed-integration/manifest.json)
retains five successful fitted/query/fixed/adaptive/donor/batch fixture checks,
six expected numerical-corruption rejections, and all six complete Kang/Hagai
checks (24,673/13,863 cells; seeds 7, 19, 41). The complete cohorts cross row-batch
boundaries. Maximum corrected-coordinate discrepancy was
2.3803181647963356e-13; maximum final-objective discrepancy was
1.5660361896152608e-11. The known biological preservation failures are unchanged.

Checker SHA256:
`cea61b0225de477e1826885c26080753466f6e6472a2acffc58a0cab4b40cc2e`.
Archive: 129 members, 134,051 compressed bytes; manifest SHA256
`e18e13793b7d049fa655bca858760bd3b91ececde84b92d025310fb9cd2e16f3`.
Verify it with the shared `../Benchmarks/HIRISA/verify_archive.py`.
The full HIRISA check is predeclared and pending, not passed by this archive.
