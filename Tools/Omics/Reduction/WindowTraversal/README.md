# Sparse PCA window traversal experiment

Full Kang H5AD PCA (24,673 cells), physical M4 Pro, 2026-09-12.
The published owner at `1fc2257098dd2f54ebaf9d10d229b5d8f19ec17e` is the baseline.
`candidate.patch` moves record mapping checks outside the inner loop while retaining
16 MiB windows, source order, FP64 accumulation, and cancellation at window boundaries.
The candidate built successfully but is **not promoted**: no useful end-to-end
speedup was established.

| Native PCA publication | Three elapsed times (s) | Median (s) |
| --- | --- | --- |
| Original | 4.030601, 4.040216, 4.005903 | 4.030601 |
| Candidate | 4.367744, 3.996794, 4.018462 | 4.018462 |

Runs used the fixed order original/candidate/candidate/original/original/candidate,
one CPU worker environment, the same source and plan, and separate fresh outputs.
Timing includes the complete CLI process but excludes subsequent comparison hashes.
The 0.3% median difference is not evidence of useful acceleration; the slower first
candidate is retained. This is six warm shared-host observations, without uncertainty
estimation, memory measurement, or a scverse PCA comparison.
All six commands succeeded. Scores, loadings and metadata match the previously
qualified full Kang output byte-for-byte. Boundary controls cover two full traversals
across a 16 MiB boundary plus three records, prior random access to the final record,
and cancellation before traversal. They do not constitute full package validation.

The prior unmodified full CLI profile is retained in `profile.txt` and
`profile-receipt.json`. Repeated project/transpose record reads appear prominently
in the sampled call tree. Inclusive sample counts overlap and include waiting
threads; they are not elapsed stage timings. The profiled 4.142 seconds includes
sampling overhead and is excluded from the comparison above.

`check.py`, `Controls.swift`, build and control logs, binary hashes and result checks
are retained. Runtime roots are `/Users/n/numivivo-pca-window-20260912` and
`/Users/n/numivivo-pca-profile-20260912`; large executables and H5AD outputs remain
there. The reproduction driver records the exact source, plan and output paths.
Apply the candidate patch to the baseline and build using
`bash Tools/Omics/H5AD/build.sh BUILD --with-cli` before running the driver.

Next: benchmark sparse project/transpose arithmetic with a bounded native kernel
or a carefully specified GPU representation. Preserve CSC as well as CSR traversal
semantics and compare complete PCA outputs before downstream qualification.
No biological prediction, million-cell PCA, or GPU PCA qualification is added here.
