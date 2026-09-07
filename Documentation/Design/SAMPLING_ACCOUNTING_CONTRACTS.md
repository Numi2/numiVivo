# Sampling accounting contracts

The native platform gate builds the complete product and test targets, executes
actual Metal fixtures, and checks the real workflow CLI. Portable numerical
tests do not replace that gate.

## Alchemical microstates (v2)

The request, result and proton-reservoir calibration bind temperature and
reference pH. `thermodynamicsRequest` may reweight to a different target pH, but
must not change the temperature or reference-pH context of the computed free
energies. A temperature change requires new free-energy/calibration evidence.
Old v1 requests/results are not silently assigned a missing calibration context.

Conditional uncertainty combines MBAR bootstrap and reservoir standard deviations
only when every required contribution is known. Unknown is `nil`, never zero.
For equal proton stoichiometry the reservoir cancels and its uncertainty is not
required. The reported root-sum-square is conditional on the documented
independent-component approximation, not a joint predictive uncertainty model.
Assumed sampling or required reservoir evidence remains assumed in the
thermodynamics adapter. Result reconstruction against the retained original
request is still required when consuming external or cached evidence.

## String refinement (v2)

Convergence requires small perpendicular force **with** a covariance-agnostic
standard-error guard and small actual final node displacement. Smoothing and
arc-length reparameterization are included in that displacement. A zero observed
force with large error bars must not pass. The component standard errors bound
RMS projected error by the triangle inequality; the configured multiplier does
not turn this bound into a confidence interval.

This updater uses the declared dimension-scaled Euclidean CV metric. It does not
infer a molecular mobility tensor or prove that the resulting path is a physical
minimum-free-energy path. Finite-restraint mean forces differentiate the
restrained center free energy; restraint-strength sensitivity remains necessary.
Equal distances along the provisional polyline do not imply equal Euclidean
chords across its corners. Regression tests use the actual arc-length
construction rather than the false chord-equality condition.

## Runtime and registry consistency

Disabling the neighbor list permits zero neighbor skin. Enabling it requires a
strictly larger radius in device precision. This does not weaken geometric
cutoff checks or neighbor overflow handling. Exact FP32 checkpoint identity
remains enforced.

The QM/MM registry module owns its free-energy adapters. The platform root must
not register them a second time or silently deduplicate conflicting definitions.
The original target-reference and finite-drug workflow interfaces remain
available alongside the newer sampling operations.
