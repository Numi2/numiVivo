# Native adjusted QL cohort hypothesis test

Frozen before native inference results. Base 9da4ee0d50aea97d6b8312c3acdd4796225c124a.
Preserve the full 12-item single-cell goal. This stage connects native abundance,
global scale/refit, adjusted residuals and robust unequal-DF moderation to a
one-constraint cohort QL hypothesis test with indexed failures and BH results.

Use all 58 original Kang, Hagai and Crowell full-support arms and all 483,576
gene/arm rows. Bind exact counts/design/offsets/contrast, original native
refitted fits, native moderation and source identities. Reuse original fitting
evidence for comparisons; test the complete new native entry point by computing
its upstream stages and requiring exact equality to the published fit inputs.
Do not substitute reference abundance, posterior scales or test probabilities.

For c'beta=0, retain the existing largest-contrast-pivot null design and fixed
scaled trend dispersion from the global refit. Compute the one-DF NB likelihood
ratio using cancelled count-only terms. The QL statistic is LR/posteriorVariance.
Denominator DF is min(adjustedResidualDF + genePriorDF, sum ordinaryResidualDF
over the complete fitted family). Preserve zero statistics and logarithmic tails,
and distinguish probability underflow from numerical failure. The Poisson bound
is not applied to this modern adjusted-residual method: pinned edgeR 4.10.5
glmQLFTest explicitly disables it; its documented bound belongs to legacy QL.
Legacy QL and active-donor borrowing are separate, still-unqualified contracts.

Qualify native constrained means against tightly converged pinned edgeR
mglmLevenberg fits to the same null design and dispersion, with relative error
at most 2e-6, dividing by max(1,abs(reference)). Independently recompute null
scores (native limit 1e-7 with 1e-10 rounding allowance), constraint residuals,
LR (2e-6 relative limit), F and denominator DF. Same-native-statistic R F log-tail
error must be at most 2e-7 absolute; posterior/F/DF/BH arithmetic error at most
2e-12 relative. Retain comparisons to the original default-reference QL study
descriptively because upstream fits and default prior optimization differ.
Report all failures and family membership; never drop a failing arm or relax a
limit to obtain a pass. Record sham calls at the original threshold without
claiming FDR calibration, truth or a superiority result.

The public cohort option must enter through the existing pseudobulk expression
owner, preserve paired/batch design and offsets, retain gene identity/statuses,
and expose QL statistics without presenting conditional Wald intervals as QL
intervals. Unsupported option combinations must fail explicitly. Keep existing
Wald/LRT behavior unchanged. Controlled fixtures qualify contracts and failure
handling only; real-family comparisons own experimental-data evidence.

Run native CPU FP64 on the physical Mac mini after workload checks. Persist
outputs remotely, transfer compressed files over IPv4, and checkpoint exact
input/executable/checker/result hashes. Retain failed attempts and typed partial
results. This stage does not promote defaults or complete biological calibration,
effect-interval inversion, GPU acceleration or million-cell qualification.
