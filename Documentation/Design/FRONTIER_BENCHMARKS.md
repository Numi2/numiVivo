# Frontier benchmark campaign v1

This is the first executable milestone of the frontier roadmap: independent
Hamiltonian comparisons and bounded realistic-system execution. It is not a
certificate for all six roadmap milestones.

## Preregistered panel and acceptance

Inputs come from pinned openmmtools revision
`f6ef22a8b9f66e582df2ffa62f3bb6516de43536` and OpenMM 8.6.0 package data.
The panel comprises TIP3P water, NaCl in water, solvated alanine dipeptide,
DHFR, the supplied vacuum T4 lysozyme/ligand complex, solvated DNA, and POPC.
The original geometry is preserved. Preparation exceptions are retained.
The original water source failed the reference engine's fixed-cutoff cell
admission. A separately identified `water-orthogonal` package fixture is added;
the original failed case stays in the scorecard and raw evidence.
This panel exercises prepared parameter input; it does not validate automatic
chemical perception, protonation or parameter fitting.

Three independent FP64 OpenMM Reference evaluations use exactly FP32-rounded
coordinates and cells: source, common translation, and deterministic perturbation.
Native and reference calculations omit LJ dispersion corrections and switching.
Periodic comparisons use a 0.8 nm cutoff; native PME spacing is 0.06 nm and
both requested electrostatic tolerances are 1e-7. Their mesh algorithms differ;
the requested tolerances do not constitute measured error bounds.

Limits fixed before the first native measurements:

- Absolute energy error per particle <= 0.002 kJ/mol.
- Force component RMS error / reference component RMS <= 0.001.
- Maximum force component error / reference component RMS <= 0.01.
- Reference RMS normalization has a 1 kJ/mol/nm floor.
- Evaluated native coordinates/cell equal the requested reference geometry.
- Accepted native checkpoint is unchanged by all static probes.

All three comparisons must pass. Optional 100-step 1 fs trajectories at 300 K
must commit every step to pass the execution smoke gate. Such trajectories
provide neither equilibrium evidence nor a throughput benchmark. Existing
analytic PME/thermal/NPT tests remain separate gates. Failed limits are never
enlarged to make a measured candidate pass.

## Interface and reproducibility

`VivoMDBenchmarkRequest` binds the complete classical system, configuration,
reference observations, declared provenance, limits and bounded dynamics work.
`numivivo md-benchmark request.json --output new-report.json` runs the owning
Metal runtime and publishes a no-clobber report. Unsupported preflight returns
an explicit report and nonzero status. Malformed input or execution errors
return nonzero with stderr; campaign orchestration must retain those failures.

Passing means agreement with the supplied observations. Descriptive provenance
inside a JSON request is not authentication of an independent calculation.
The campaign must hash its source inputs, serialized OpenMM System, generator,
request, report, executable, shaders and logs. Native reports additionally bind
the canonical request fingerprint and numerical execution contract.

The optional `dynamicsPreparation: "projectConstraints"` invokes the owning
transactional position/velocity projection before assigning thermal velocities.
Absent/`preserve` retains prior behavior. The prepared checkpoint is explicit
in the dynamics result; static reference coordinates are never projected.

Compensated positions introduced by the [v5 position profile](MD_NUMERICAL_PROFILE_V5.md)
remain opt-in under the [v6 RATTLE profile](MD_NUMERICAL_PROFILE_V6.md), through
`configuration.positionPrecision: "compensated"`. The independent Python
verifier checks the saved position/velocity constraints and exact checkpoint
words. Optional `dynamicsObserveEvery` records at most 1002 energy/temperature
observations at accepted steps. The initial matched-duration NVE study fixes
its limits in `Tools/Benchmarks/nve_policy.json` before measurement.

The separately named smooth-LJ panel recalculates independent observations
after enabling a 0.7 to 0.8 nm switch in both engines. The original sharp-cutoff
model and preparation failures remain retained. Its measured conservation and
refinement results are in the [frontier audit](../Audit/2026-09-08_FRONTIER_MD_REFERENCE_PANEL.md).

The reference environment is isolated and never used as a production force
provider. `Tools/Benchmarks/requirements.txt` pins its packages. Source/package
data are downloaded into the external evidence directory, not silently added
as unlicensed copies to the repository. Upstream distribution terms remain
applicable; links and identities accompany generated observations.

## Subsequent roadmap gates

The [next ensemble milestone](MD_ENSEMBLE_QUALIFICATION.md) now provides a fixed
5 ps conservation campaign and independent kinetic/configurational NVT analysis.
Its [measured follow-up](../Audit/2026-09-08_INTERACTING_ENSEMBLES.md) passes all
seven prepared systems at three time steps over 5 ps, including independent
endpoint energy checks. The completed 50 ps NVT matrix passes kinetic checks at
0.5 fs, fails kinetic checks at 1 fs, and lacks sufficient configurational sampling
at both steps. A separately declared 210 ps finer-step follow-up is running.

1. Extend this panel to interacting NVE/NVT/NPT ensemble distributions,
   timestep refinement, constrained motion and numerical failure recovery.
2. Qualify reduced triclinic execution across all geometry-dependent owners.
3. Measure bounded multi-step execution, neighbor reuse, constraints and PME;
   compare matched-accuracy simulations on the same Apple hardware.
4. Extend existing umbrella/MBAR owners with restartable replica exchange and
   neutral FreeSolv hydration campaigns.
5. Qualify existing factorized chemistry/analytical derivatives and evaluate
   learned forces against held-out static and trajectory evidence.
6. Publish realistic end-to-end workflows and a release qualification panel.

References: [physical validation](https://doi.org/10.1371/journal.pone.0202764),
[OpenMM Metal](https://github.com/philipturner/openmm-metal),
[GROMACS performance](https://manual.gromacs.org/current/user-guide/mdrun-performance.html),
[FreeSolv](https://github.com/MobleyLab/FreeSolv),
[2026 learned-potential comparison](https://arxiv.org/abs/2601.16331).
