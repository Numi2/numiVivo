# Applying the declared reference tolerance to the owning solver

The first sensitivity attempt failed its method constraint: pinned
`glmFit.default` does not forward `tol` or `maxit` from `...` to its likelihood
fit. It calls `mglmLevenberg` with maxit=250 and the solver's default tolerance.
All ten attempted outputs and the adapter source remain under
`ignored-tolerance-attempt`; their unchanged discrepancies are not a completed
tighter-optimization experiment.

The repaired adapter calls the same exported `mglmLevenberg` owner directly
for full fits already using that solver, with tol=1e-10, maxit=1000 and the
default fit as its starting coefficients. Replace only returned fitted-model
fields in the existing DGEGLM object before ordinary glmLRT. Preserve its
counts, offsets, dispersions and design. One-way full fits and glmLRT's null
fits retain their existing paths. Record solver identity and whether tolerance
was actually refined for every gene. Repeat all ten cases and preserve both
default-reference failures and the failed first sensitivity attempt.
