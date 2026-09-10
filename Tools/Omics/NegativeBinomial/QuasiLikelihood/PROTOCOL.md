# Quasi-likelihood stages before native cohort integration

Freeze before computing the new stage results. Reuse all twenty original
Kang/Hagai default-support sham analyses and nine full-support treatment
analyses (Kang default, Hagai default, seven Crowell populations) from the
published LRT source manifest. Exclude active-donor policy analyses explicitly:
QL prior borrowing across heterogeneous retained designs needs a separate
native contract. Do not describe these 29 analyses as full support-policy
coverage or untouched validation. Keep every original tested gene, pseudobulk,
design, contrast, count and offset; no outcome-dependent selection.

Run four pinned edgeR 4.10.5 reference arms for every case: native stored
abundance-trend dispersion versus edgeR robust estimated trend on the same
tested family; each with modern adjusted-deviance QL (legacy=FALSE) and legacy
QL (legacy=TRUE). Use robust=TRUE, abundance.trend=TRUE, prior.count=0 and
winsor.tail.p=(0.05,0.1). NB dispersion is a cohort trend, not native per-gene
MAP dispersion. Keep native offsets in effective-library units by adding the
same geometric library-count constant to every sample. This only changes the
intercept; use that identical offset in all four arms. Do not run TMM here.
Record all counts/feature IDs, requested and actual dispersions, QL average
scale, raw/adjusted deviances and residual DF, unit adjustments/leverages,
prior DF/scale, posterior scale, LR, F, denominator DF, p-values, BH, optimizer
flags/scores, warnings, messages and errors. Record dispersion caps as failed
exact-dispersion method constraints while retaining outputs. Never tune the
fit or filtering to reduce sham calls.

Independent checks reconstruct F=LR/posteriorScale and denominator DF as
min(priorDF+residualDF,sum(n-p)), the F survival probability and BH. Legacy
Poisson bounds are retained explicitly and reproduced as max(F-tail,bound).
Tolerance: F/DF relative 1e-10, p/BH absolute 2e-10; posterior arithmetic
relative 1e-10 when prior DF is finite. Modern unit adjusted DF/deviance sums
must agree with aggregate fields to relative 1e-9. Preserve every failed check,
nonfinite value and unavailable gene. Fit stopping diagnostics are observations,
not proof of strict convergence from a false package failure flag.

Add a native, stable NB unit-deviance primitive needed by residual QL modeling.
Its count domain is the existing exact UInt64-through-2^53 domain; support
zero means only for zero observed count, the exact Poisson limit, and finite
NB dispersions through 100. Check Poisson/zero/near-saturated/extreme cases
against high precision and all actual reference fitted-count pairs against
independent arithmetic (absolute 2e-7 plus relative 2e-10). This arithmetic
primitive does not itself implement QL moderation or calibrated inference.

Report all four complete families, original native Wald/LRT call counts and
stage changes. Overlapping sham splits are not independent experiments and
treatment discoveries are not known power. These reference-stage results guide
native adjusted-deviance and robust prior ownership; neither a native QL method,
FDR-control claim nor a default promotion follows from this study. Retain
source/script/runtime hashes and externally stored complete matrices. Avoid
dense cells-by-genes arrays; only the small donor pseudobulk matrix is dense.
