# Numerical repair before promotion

After the initial complete cohorts and their biological diagnostics passed,
an extended near-underflow check found a 4.5777e-5 delta error at exponent 740.
Total weights agreed but BLAS subnormal product accumulation differed from scalar
arithmetic. Preserve that failed check and the initial source/executable/results.

For query total weight below 1e-280, recompute the complete Gaussian average with
original-order direct scalar distance, exp, numerator and denominator sums. Charge
every recomputed distance term to maximumWork; check cancellation every 256 anchors.
No rescaling, cutoff, dropped anchor, biological parameter or margin changes.
This is a numerical defect repair, not tuning to the biological scores.

Run the expanded 31 numerical cases, admission/error cases and sanitizers. Rebuild
and rerun the same three complete cohorts with scalar plus tiled plus replay.
If repaired full scores/anchors are byte-identical to the frozen first trial,
reuse its complete biological evaluation with explicit payload identity checks.
Repeat the full Hagai publication/verification path. Preserve legacy scalar plans.
