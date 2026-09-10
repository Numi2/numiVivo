# Adaptive discrete deviance-moment summation

Frozen before native adaptive results. Base: 7c713b8f8a84d09d2dc41bf86cfcd198220c11c6.
The direct mode-centered engine and all its work-limit failures remain retained.
The goal is a practical native QL moment stage for the complete existing family,
not a change to fitted means, dispersions, feature support or statistical outcomes.

Add an explicitly selected adaptive summation method. Preserve direct summation
as the existing default. The adaptive selection uses direct summation when mean
is below 10,000 or dispersion is zero; otherwise it evaluates integer endpoints
of count blocks. It uses cubic Hermite interpolation of f(k)=P(k), P(k)D(k), and
P(k)D(k)^2. For an integer block [a,b], n=b-a, its polynomial sum is
(n+1)(f(a)+f(b))/2 + (n^2-1)(f'(a)-f'(b))/12.
The absolute remainder is bounded by sup|f''''| (n^5-n)/720. Two-point and
one-point blocks are exact. Subdivide until explicit mass/first/second moment
error budgets pass; work exhaustion is an error, without an asymptotic fallback.

Bound f'''' using the product rule, exponential derivative polynomials for the
NB mass, monotonic log-PMF slopes, polygamma-series derivative bounds, and exact
monotonic bounds on deviance and its derivatives. For n>=2, with z=a+min(1,r),
|log(P)^(n)| <= |r-1| ((n-1)!/z^n + n!/z^(n+1)). Retain geometric left/right
bounds on omitted mass and deviance moments. Propagate separate summation and
tail error into normalized mean and variance, and require their combined bounds
to meet the original relative tolerance 1e-10. Bounds concern truncation and
interpolation; floating-point rounding is not interval-certified.

Evaluate the PMF with stable saturated likelihood and Stirling correction;
never subtract large log-gamma values. Use stable digamma differences for the
continuous derivative at integer endpoints. This is finite discrete summation,
not replacing the NB distribution with a continuous Gamma distribution. The
Hermite interpolation remainder follows the standard repeated-node polynomial
remainder; the needed gamma/psi identities are in NIST DLMF sections 5.11/5.15.

Qualification:
- Retain the original 77-point grid and focused native tests. Add polynomial
  sum/remainder tests and independent checks of PMF, slope and derivative bounds,
  budget failures and method dispatch. Add an extended grid with means
  {10000,100000,1000000} and dispersions {0.001,0.03,0.1}, plus mean 10000 with
  dispersions {0.7,1,4}. Compare with independent direct PMF summation in bounded
  chunks at relative tolerance 2e-7 for mean, variance, A and K. This tolerance
  covers ordinary floating arithmetic separately from the requested analytic
  approximation bound.
- Rerun all twenty originally unavailable gene/arm residual inputs with the
  adaptive method, default one-million evaluation limit, unchanged tolerance
  and all original counts/means/designs/scales. Compare all 118 existing higher-
  work direct moments, and independently check all 124 adaptive moments. Retain
  any remaining failure. Measure process time/RSS and evaluation work.
- Apply the explicit adaptive method to both modern QL arms for all 29 existing
  cases. Compare all previously available direct results, and check every newly
  available moment independently. Record the selected summation method and
  bounds. No smaller family can establish complete coverage.
- Measure the pilot before scheduling the entire family. Do not infer global
  QL scale estimation, robust moderation, QL hypothesis tests, FDR calibration,
  general performance superiority or million-cell execution from this stage.
