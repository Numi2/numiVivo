# Classical MD numerical profile v4

`numivivo.org/md-metal-numerics/v4` corrects the sign of periodic-torsion
forces. The torsion energy is unchanged. The previous shader contracted
`dE/dphi` with coordinate factors requiring `-dE/dphi`, reversing all four
atomic force contributions. Isolated OpenMM Reference evaluations at phases
0, 0.37 and pi agree with independent central energy differences; v3 fails
the force comparisons while its energies agree.

This changes the dynamics and minimization of systems containing nonzero
torsional forces. Historical v3 trajectories must not be relabeled as v4.
The existing checkpoint, sampling and QM/MM identity owners inherit the new
contract. Old checkpoints remain historical data; direct restore rejects.
Explicit import of physical state starts a new numerical run with retained
source provenance. The algorithm change does not retroactively validate old
torsional dynamics.

`VivoMDMetalRuntime.projectConstraints()` is a separate explicit preparation
operation. It projects candidate positions and velocities using the existing
constraint kernels, checks the candidate, and publishes one checkpoint without
advancing the MD clock. Failed projection or cancellation preserves the last
accepted state. It is neither energy minimization nor equilibration.

The frontier benchmark request can opt into this preparation for dynamics.
Static energy/force comparisons always retain the original reference geometry.
The prepared thermal checkpoint and final accepted checkpoint distinguish
preparation from subsequent trajectory execution.

Required evidence includes independent isolated torsion forces and derivatives,
the realistic reference panel under unchanged limits, transactional preparation
and old-contract rejection, existing MD/PME/thermal/NPT regressions, and the
full native regression suite. Performance and equilibrium claims remain
separate from static agreement and short trajectory completion.
