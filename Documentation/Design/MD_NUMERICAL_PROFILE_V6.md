# Classical MD numerical profile v6

The preregistered realistic NVE study rejected v5 on all seven prepared systems.
Energy decreased approximately linearly with time, and halving the time step
only reduced the error by about a factor of two. Satisfying distance and velocity
constraints alone did not establish correct constrained dynamics.

Profile v6 corrects the constrained drift to use RATTLE position impulses along
the constraint normals at the **start of the drift**. Each nonlinear iteration
uses the current distance residual and the current/reference distance dot
product in its denominator. The original closest-point projection used the
changing trial normal, damping tangential motion. Velocity corrections retain
the full constrained displacement divided by the drift duration, followed by
the final tangent velocity projection.

Explicit initial projection and minimization retain their geometric projection;
they are not dynamics steps. NVE and both Langevin-middle half drifts use the
fixed-reference RATTLE path. The change applies to FP32 and compensated modes.
All coordinates, velocities, correction words and clocks retain the existing
transaction boundary. Nonpositive correction denominators reject the candidate.

A reusable start-of-drift reference costs 16 bytes per particle in both modes.
Compensated mode now has six correction buffers (96 additional bytes per
particle), checked against the device buffer and remaining working-set limits.
Force precision, cell restrictions and unsupported compensated paths retain
the v5 contract. Old numerical profiles cannot directly resume under v6.

Required evidence includes a free rigid rotor's energy/angular momentum, the
same realistic reference and constrained-execution panel, and the unchanged
three-timestep NVE policy. The failed v5 study is retained as the baseline.
Passing these short studies would not establish interacting equilibrium,
transport properties, general bitwise PME replay or performance leadership.

Reference: Andersen, [RATTLE: A velocity version of the SHAKE algorithm](https://doi.org/10.1016/0021-9991(83)90014-1),
Journal of Computational Physics 52 (1983), 24–34; maintained
[LAMMPS SHAKE/RATTLE documentation](https://docs.lammps.org/fix_shake.html).

The original periodic panel uses a sharp LJ cutoff, so cutoff crossings can
produce energy jumps independent of the integrator. The initial v6 study retains
that model and its failed/inconclusive refinement results. A separately named
smooth-LJ panel switches interactions from 0.7 to 0.8 nm in both native and
OpenMM models and recalculates all independent observations. It changes the
Hamiltonian; it does not relabel the original panel or relax the NVE policy.
See [OpenMM's switching-function definition](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html).
