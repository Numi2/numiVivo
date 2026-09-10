# Native abundance covariates and QL global-scale integration

Freeze before new native/reference results. Base d45f93bb2b84e51208a4dcc9dc0c8cd7f78fe6db.
This stage replaces externally supplied abundance coordinates in an explicit
global QL entry point. It does not change cohort Wald/LRT defaults or add QL
priors, posterior moderation or hypothesis tests.

For counts y_i, supplied log library offsets o_i, scalar dispersion phi and
average prior count p (default 2), define L_i=exp(o_i), a_i=p*L_i/mean(L),
z_i=y_i+a_i and L'_i=L_i+2*a_i. Fit an intercept-only NB2 model with means
L'_i*exp(beta). The score is sum((z_i-mu_i)/(1+phi*mu_i)); its derivative is
strictly negative for positive total augmented counts. Report
(beta+log(1e6))/log(2). Phi=0 is the Poisson limit. Positive-prior zero-count rows
have finite abundance; zero total augmented counts have no finite log abundance.

Implement from this mathematical definition and primary documentation. Use
centered offsets, stable prior scaling, a sign bracket and bounded root search.
Preserve score, bracket-width and evaluation diagnostics. Reject invalid inputs,
unrepresentable scaled priors and exhausted work explicitly. Raw count entries
must remain exact UInt64 values through input validation; no pseudo-counts are
added to the original NB family response.

Compare controlled cases covering zero counts, prior counts 0/0.5/2/10, Poisson
and NB dispersions, unequal libraries, offset shifts, permutations, large exact
counts, and invalid/work-limit cases. Use pinned edgeR 4.10.5 aveLogCPM and an
independent scalar score root in Python; retain disagreements. Native abundance
agreement is required within 2e-8 absolute log2 CPM on the supported comparison
domain. For challenging controlled cases, independently check with high
precision arithmetic and preserve reference convergence limitations.

For both original full-support trend arms of all 29 Kang, Hagai and Crowell
cases, first compare abundance with the already retained d45f93bb stage inputs.
Then run the explicit native-abundance global-scale path on unchanged counts,
designs, offsets, contrasts and trend dispersions. Require all 58 arms, no
silent family exclusion, and complete structured results. Compare resulting
global scales, refit means and final adjusted residuals with the retained
supplied-abundance native stage within 2e-7 relative to max(1,abs(reference)).
If a smoother branch amplifies a small abundance difference beyond that limit,
retain and independently arbitrate it rather than weakening the limit.

Check native abundance score roots independently. Native final GLM scores must
still satisfy 1e-7, with only 1e-10 absolute independent-rounding allowance.
Retain all failure receipts, source/runtime hashes and unchanged prior input
hashes. Measure a pilot, confirm host workload ownership, checkpoint completed
outputs, and keep large matrices external by exact hashes. Native computation
is CPU FP64 on the physical Mac mini. These numerical and pipeline checks do
not establish statistical calibration, biological truth, general throughput
superiority, GPU acceleration, or million-cell execution.
