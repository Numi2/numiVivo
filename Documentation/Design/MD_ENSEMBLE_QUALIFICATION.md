# Interacting NVT qualification and longer NVE checks

This extends the v6 force/constraint milestone using the same owning Metal
runtime. Production force evaluation stays native. Independent Python analysis
reads retained observations; it does not provide forces to the simulation.

## Fixed campaigns

The separately named smooth-LJ panel retains the previous Hamiltonian and static
limits. An explicit 32-iteration execution variant uses the runtime's existing
default, with the same 1e-6 constraint tolerance. Short comparisons and NVE
refinement check this variant before longer runs. A failed variant remains
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

References: [Merz and Shirts, physical validation](https://doi.org/10.1371/journal.pone.0202764),
[OpenMM's documented velocity convention](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LangevinMiddleIntegrator.html).
