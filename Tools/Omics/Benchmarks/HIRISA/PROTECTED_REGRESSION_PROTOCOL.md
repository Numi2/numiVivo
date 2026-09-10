# Fixed-membership protected-stratum regression experiment

Declared 2026-09-11 after the published HIRISA annotation-retention failures.
This is a single development hypothesis on inspected data, not independent
validation. No native integration algorithm or default changes in this experiment.

Use every original cell and all twenty PCs. Retain the exact native seed-7 final
soft cluster memberships, donor identities and original PCA. Reconstruct the
original fixed-ridge donor correction first. Then fit exactly one alternative:
within each cluster, an unpenalized intercept for each original experimental
preparation × treatment stratum, plus donor effects with the unchanged ridge
penalty one. Subtract only donor effects, weighted by the same memberships.
Use the original active-donor criterion and all original cells. Do not use author
cell-type labels, annotation confidence, RNA program scores, held-out response
outcomes or the new diagnostic scores to fit this correction.

For cluster c, let m_sb be membership mass and t_sb the membership-weighted PCA
sum in stratum s and donor b. Eliminate stratum intercepts to solve:

    G = diag(sum_s m_sb) - sum_s outer(m_s,m_s)/sum_b m_sb + I
    rhs_b = sum_s t_sb - sum_s m_sb * sum_b t_sb / sum_b m_sb
    G beta = rhs

Each stratum intercept is its weighted PCA mean minus the mass-weighted donor
effects. A stratum represented by only one active donor supplies no contrast to
identify that donor effect. Ridge keeps the effect system nonsingular; it does
not turn confounded data into evidence. Check every cluster against independent
augmented least squares over weighted stratum/donor means, and check normal
residuals and complete corrected coordinates.

Freeze source hashes, stratum mapping, protocol and implementation before fitting.
Freeze both reconstructed and protected coordinate files before running metrics.
Apply the unchanged full-cohort annotation protocol, including all 114 folds,
471 comparisons, insufficient strata and original mean/per-fold margins. Use the
existing response and within-library program diagnostics unchanged as separate
endpoints. Report conditional donor-associated variance and mean-drift descriptions
so retention cannot be mistaken for useful correction merely because effects
shrink. Do not select correction strength or margins from these results.

The baseline and alternative share frozen native memberships. This isolates the
regression step; it is not a refit of soft clustering, a new native execution,
prospective mapping, proof of batch-effect identifiability or biological truth.
A favorable result would justify a separately implemented and qualified iterative
native method. A failure remains a failure; no parameter sweep is allowed here.
