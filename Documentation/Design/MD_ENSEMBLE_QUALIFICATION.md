# Interacting NVT qualification and longer NVE checks

This extends the v6 force/constraint milestone using the same owning Metal
runtime. Production force evaluation stays native. Independent Python analysis
reads retained observations; it does not provide forces to the simulation.

## Fixed campaigns

The separately named smooth-LJ panel retains the previous Hamiltonian and static
limits. An explicit 32-iteration execution variant used the runtime's existing
default, with the same 1e-6 constraint tolerance. It failed admission on several
systems and remains failed. A 64-iteration variant passed all seven short force,
constraint and NVE-refinement cases; the prospective NVT policy uses 64 iterations.
This increases solver work without changing its accuracy requirement. A failed variant remains
failed; solver or scientific limits are not loosened to promote it.

- `long_nve_policy.json` extends conservation/refinement from 0.1 to 5 ps, using
  the same energy/error policy and three time steps. This is 50 times longer,
  still far shorter than a production molecular research trajectory.
- `ensemble_policy.json` fixes rigid orthogonal water, 300/305 K, three independent
  seeds, two time steps (1/0.5 fs), friction 10/ps and 50 ps per run. The first
  10 ps is excluded prospectively. Observations every 0.05 ps retain 800
  production samples per seed; one-ps blocks support correlated uncertainty
  estimates. This elevated friction is a declared test condition, not a transport
  model or a claim about every thermostat setting.

The higher-temperature seeds add the fixed offset 1,000,000 to the three listed
base seeds. Thus the two temperature samples use disjoint random streams, while
each temperature/seed pair retains the same prepared state across time steps.
This was fixed before any NVT measurements; early policy snapshots remain visible.

The subset and excluded prepared cases are explicit in derived manifests.
Original preparation failures remain in every derived campaign. NVT input here
is rigid three-site water with independent triangular constraints and no fixed
atoms or center-of-mass removal. Its kinetic degrees of freedom are independently
counted as six per molecule. Other molecular topologies require their own count.

## Independent analysis

At equilibrium, total kinetic energy has a gamma distribution with shape
`degreesOfFreedom / 2` and scale `kB * targetTemperature`. The count and target
come from the input model, not a fit to observed data. Check mean, width and
distribution shape, retaining block-bootstrap uncertainty and sample correlation.
Constant traces, insufficient effective samples, unresolved block correlation,
replica disagreement or drift cannot certify an ensemble.

The ratio of potential-energy distributions at two temperatures has an expected
logarithmic slope `1/(kB*Tlow) - 1/(kB*Thigh)`. A free-intercept logistic fit checks
this slope without needing an analytic configurational density. Independently
resampled blocks preserve temporal correlation within each replica and temperature.
Insufficient overlap or precision is inconclusive. Moment, distribution and slope
rejections remain failed observations. Every prescribed time step must be assessed.

Native v6 observations use on-step BAOAB velocities. OpenMM's LFMiddle reports
half-step velocities, so its raw kinetic energies must not be treated as the same
estimator. Finite-step kinetic bias and the finer-step result remain visible.

These are necessary kinetic/configurational consistency checks. Passing does not
prove complete equilibration, every degree-of-freedom partition, experimental
water properties, NPT, free energies, transport or general scientific leadership.
Policies, source/binary/shader identity, exact commands, all inputs, reports and
failures must accompany measured results.

## Declared sampling extension

The first completed 300 K, 1 fs subset had only 142.4 total effective potential
samples, below 200, and one-ps block correlation 0.419, above 0.3. Its kinetic
mean also missed the fixed 3.5-sigma threshold (3.579); that observation stays
failed. These are interim subset findings, not a completed matrix verdict.

`ensemble_continuation_policy.json` declares a separate adaptive follow-up before
any continuation data: all six original 0.5 fs replicas continue to 210 ps. The
first 10 ps remains excluded; five-ps blocks provide 40 complete blocks per
replica. Every statistical acceptance limit remains inherited unchanged. The
longer policy is explicitly informed by the initial coarse subset, not a blinded
preregistration or a retroactive relabeling of the 50 ps result.

`continue_ensemble.py` requires the original complete matrix and qualification,
retaining failed and inconclusive outcomes. It checks each selected parent's
passed native/static execution and exact input hashes, and restores its existing
final checkpoint with the identical executable and complete shader bundle.
A zero-step restore must preserve all checkpoint words, velocities, cell,
fingerprints, time and stochastic step counter. Ten-ps chunks preserve the same
configuration; no velocity initialization or new equilibration exclusion occurs.

The existing `md-run` command emits observables together with sampled JSON
positions. Both sampling intervals must be supplied. JSON snapshots reconstruct
the compensated high/low positions in Double; no FP32 trajectory archive is used.
Every chunk retains commands, raw reports, exact checkpoints, hashes and an
independent final kinetic/potential audit. Offline analysis verifies the complete
chain, rejects missing/duplicate observations and assesses all original plus new
production observations with the declared wider blocks. This does not claim
bitwise PME trajectory replay or extend the result to other time steps.

The [completed measurement audit](../Audit/2026-09-08_INTERACTING_ENSEMBLES.md)
records the original 50 ps failures and inconclusive results and the subsequent
prepared-case pass at 0.5 fs over 210 ps. The original failed water preparation
remains retained; the broader frontier roadmap stays open.

References: [Merz and Shirts, physical validation](https://doi.org/10.1371/journal.pone.0202764),
[OpenMM's documented velocity convention](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LangevinMiddleIntegrator.html).
