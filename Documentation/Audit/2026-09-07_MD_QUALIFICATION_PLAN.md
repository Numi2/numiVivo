# MD numerical qualification plan — 7 September 2026

This is a source audit and a frozen native acceptance plan, not a record of measured numerical results. The inspected production source is `d9c9de49c069dcd0635acf6fbb43d006e4cd1ec0`; the owning classical MD algorithms are unchanged in publication `496c99aa02499832781d73192a90608c70301c7d`. No build, GPU job, trajectory, force measurement or ensemble measurement was performed for this audit. The numerical limits below were frozen before native execution; they must not be widened after observing a failure without a separately documented scientific decision.

The highest-value next gate is charged classical PME energy and force agreement against an independent direct Ewald sum. The existing execution tests do not provide that evidence. Constraint and ensemble distribution checks follow it; a passing force fixture does not qualify a thermostat or barostat.

## Supported source contract and present evidence

| Area | Executable source contract | Existing evidence and remaining gap |
| --- | --- | --- |
| Classical electrostatics | Cutoff, reaction field and cubic charge PME. PME requires neighbor-list execution and a periodic cell. | The named PME smoke fixture has zero charges on every particle. No independent native classical-PME force comparison was found. |
| Cell geometry | Orthogonal periodic cells, with minimum-image cutoff and neighbor-radius validation. | Triclinic structure storage and the multipole reference's triclinic support do not qualify triclinic MD. |
| Constraints | Iterated position and tangent-velocity projection for distance constraints between massive particles. Reported thermal degrees of freedom subtract the constraint count and assume independent constraints. | The populated native fixture contains two disjoint constraints; one-step commit does not qualify coupled constraint iteration, long-time residuals or constrained velocity distributions. |
| NVE | Verlet-style half kicks and constrained drift. | A two-particle harmonic native fixture checks energy, drift, momentum and exact restart. This is useful scoped evidence, not the full force matrix. |
| NVT | Langevin middle: half drift, stochastic velocity update, tangent projection, half drift. | The native smoke checks finite thermalized velocities and one committed step. No independent finite-time or stationary distribution assertion was found. |
| NPT | Langevin dynamics with molecular-center isotropic Monte Carlo volume proposals, symmetric in log volume. | The native smoke checks one commit and the existence of a barostat certificate. It does not check the proposal score, volume measure or equilibrium distribution. |

Owning source: [configuration and capability notes](../../Sources/NumiVivoKit/MD/VivoMDTypes.swift), especially lines 155–185; [execution preflight](../../Sources/NumiVivoKit/MD/VivoMDExecutionPreflight.swift), lines 5–47; [native step and thermalization](../../Sources/NumiVivoKit/MD/VivoMDMetalRuntime.swift), lines 180–257; [constraint and Langevin kernels](../../Sources/NumiVivoShaders/Resources/NumiVivoMDKernels.metal), lines 29–33. Velocity-rescale, Drude dynamics and unsupported polarization are blocked; explicit supported variational providers form a separate qualification scope. Capability reporting already describes the PME tolerance as a planning heuristic rather than a measured force-error bound.

The specific coverage distinction is important:

- [AppleExecutionTests.swift](../../Tests/NumiVivoIntegrationTests/AppleExecutionTests.swift), lines 32–36, assigns zero charge to every populated-fixture particle. `NVTAndPMEExecuteOnMetal` at lines 204–219 and `NPTAndMinimizationExecuteOnMetal` at lines 221–230 therefore do not exercise charged PME forces. The harmonic NVE/reference test is at lines 154–185.
- [MultipolePMETests.swift](../../Tests/NumiVivoIntegrationTests/MultipolePMETests.swift), line 16, does compare multipole mesh forces and energy with direct Ewald. It uses the sixth-order multipole assignment. Classical MD uses separate cubic spread, influence and gather code in [NumiVivoMDPME.metal](../../Sources/NumiVivoShaders/Resources/NumiVivoMDPME.metal), lines 97–148; the multipole implementation starts at line 153. Sharing FFT stages does not establish equality of these complete force paths.
- [PeriodicQMMMHamiltonianTests.swift](../../Tests/NumiVivoIntegrationTests/PeriodicQMMMHamiltonianTests.swift), line 152, exercises a charged native BO/PME path, but asserts execution, motion, coarse energy and provider/restart behavior rather than classical force agreement. [AdaptiveQMMMHamiltonianTests.swift](../../Tests/NumiVivoIntegrationTests/AdaptiveQMMMHamiltonianTests.swift) uses neutral classical particles and checks finite/repeated Hamiltonians and execution. Those tests do not close this gap.

These are qualification gaps. This audit does not establish that the corresponding runtime algorithms are numerically wrong. Prior build/reliability results remain attached to their own source and scope in the [completion reliability record](2026-09-07_COMPLETION_RELIABILITY.md).

## Priority 1: charged classical PME against direct FP64 Ewald

Use the existing public [Hamiltonian probe](../../Sources/NumiVivoKit/MD/VivoMDMetalRuntime.swift), line 388. It evaluates the same force kernels as dynamics without committing a new accepted state, time or cell. Every supplied physical coordinate must be exactly representable in FP32. Generate explicit binary-fraction coordinates, or explicitly round once before constructing both the native geometry and the oracle. Compare using the probe's returned normalized geometry and require the accepted checkpoint to remain unchanged.

### Declared fixture family

The candidate `ClassicalPMEConformanceTests.swift` uses four massive particles, charges `[1, -1, 0.5, -0.5]`, zero Lennard-Jones coefficients, no bonded forces, constraints or force provider, and a 2 nm orthogonal cell. Use a 0.75 nm real cutoff, 0.125 nm neighbor skin, PME planning tolerance `1e-7`, and fixed `16^3`, `32^3` and `64^3` meshes. Base positions in nm are `[43,73,107]/256`, `[477,101,131]/256`, `[233,85,115]/256` and `[309,351,409]/256`. They are exact FP32 coordinates away from grid nodes even on the finest grid. The nearest-image 0–1 pair crosses the box boundary at about 0.337 nm; 0–2 is near, but inside, the real cutoff at about 0.7443 nm; other pairs exercise the reciprocal-dominated region.

The exact native bound is nine cases times three grids, or 27 Hamiltonian probes: neutral base, neutral displaced, neutral translated within the cell, excluded pair 0–1, half-scaled pair 0–1, charged `[1,-0.5,0.5,-0.25]` with net charge 0.75, the charged half-scaled case, neutral whole-lattice translation, and neutral affine scaling by `17/16` into a cell with side 2.125 nm. The last case evaluates the alternate geometry/cell while retaining the accepted initial cell of 2 nm and the fixed mesh. The frozen test must record the exact displacements and translations before execution. The excluded/scaled pair remains strictly inside the real cutoff.

### Independent oracle and conventions

Use [VivoPeriodicElectrostatics.evaluate](../../Sources/NumiVivoKit/QMEnv/VivoPeriodicElectrostatics.swift), line 139, with `reciprocalOperator: nil`: a direct bounded FP64 lattice sum, not either Metal mesh implementation. Set every dipole and quadrupole to zero. The API receives particle positions in Bohr and the cell in nm; it converts the cell internally. Let `a0 = VivoAtomicUnits.bohrInNM` and `Eh = VivoAtomicUnits.hartreeInKJPerMol`.

1. Obtain the declared Ewald beta from [VivoPMEPlan.make](../../Sources/NumiVivoKit/MD/VivoPMEPlan.swift), using the native 0.75 nm cutoff and tolerance. Set oracle `alphaPerBohr = betaPerNM * a0`. Fixed meshes and reference-convergence variations must not change this split. Retain the plan's Double beta in the reference and record the native Float beta separately.
2. The primary reference uses reciprocal half widths `[18,18,18]` and `realCutoffBohr = 1.125 / a0`. Check `[14,14,14]` against `[18,18,18]` at the 1.125 nm reference cutoff, and check a 0.75 nm against a 1.125 nm real cutoff at `[18,18,18]`. Each comparison must have converted energy difference at most `1e-5` kJ/mol and maximum Cartesian force-component difference at most `1e-4` kJ/mol/nm. These checks separately bound reciprocal truncation and the real-space tail omitted by the native cutoff; 1.125 nm is an oracle cutoff, not a supported native minimum-image cutoff for this cell. If either reference check fails, the native comparison is inconclusive; revise the reference extent through a separately reviewed fixture change. Do not silently accept an unconverged reference.
3. Match the numerical Coulomb convention explicitly. MD uses `C = 138.935456 / relativeDielectric`, whereas atomic-unit conversion gives `Eh * a0`. Define `s = C / (Eh * a0)`. Convert reference energy by `Eh * s` and reference force by `(Eh / a0) * s`. A difference in constants must not be mislabeled as mesh error. Native beta, coordinates and force accumulation remain FP32 and are covered by the declared native error budget.
4. Use `.requireNeutral` for the neutral cases and `.uniformNeutralizingBackground` for the charged case. Native [PME background energy](../../Sources/NumiVivoKit/MD/VivoPMEEngine.swift), lines 75–77, is `-pi * C * Q^2 / (2 * beta^2 * V)`. The reference includes the same tin-foil uniform-background convention. Self and background terms belong in the total-energy comparison even though background adds no particle force.
5. Translate `VivoNonbondedException` charge scales into reference `primaryPairScales`. The reference applies `(scale - 1) * C * qi * qj / r` to the nearest primary pair, not to every periodic image. Native real space scales `erfc(beta*r)/r`; [its exception correction](../../Sources/NumiVivoShaders/Resources/NumiVivoMDPMECorrections.metal), line 31, adds `(scale - 1) * erf(beta*r)/r`. These sum to the same correction for a primary pair inside the cutoff.

For an exception outside the real cutoff, the native omitted real-space tail leaves a difference of magnitude `abs((scale - 1) * C * qi * qj) * erfc(beta*r)/r` relative to the reference's bare primary correction. Its energy is bounded by that prefactor times the planning tolerance when `r >= cutoff`. The force-tail magnitude is the same prefactor times `erfc(beta*r)/r^2 + 2*beta*exp(-(beta*r)^2)/(sqrt(pi)*r)`. A later outside-cutoff exception case must retain the independent bare-pair oracle and declare this truncation allowance explicitly; it must not change the oracle to imitate the kernel. Avoid exact cutoff equality because the native and reference inclusion comparisons differ there.

### Proposed predeclared acceptance

For each case, define `R` as the RMS of all `3*N` reference Cartesian force components and `D = max(R, 1 kJ/mol/nm)`. Every declared `64^3` case must satisfy `abs(E_native - E_reference) <= 0.05` kJ/mol, force-component error RMS divided by `D` at most `1e-3`, and every individual Cartesian force-component error divided by `D` at most `3e-3`. Report the denominator, raw component errors and maximum per-particle vector error. The denominator floor makes the force limit absolute for weak-force cases; it is not a guarantee of relative accuracy near zero. These limits are proposed product gates, not observed precision or a universal consequence of the input tolerance.

Grid trend is a separate, explicitly scoped test over this exact nine-case set. Compute RMS across cases of the energy errors and, separately, RMS of the per-case normalized force errors. Require aggregate `64^3` energy error to improve on `16^3` unless both are below `1e-3` kJ/mol; require aggregate `64^3` normalized force error to improve on `16^3` unless both are below `1e-5`. The per-case fine absolute gates apply regardless of the trend floor. Retain all `32^3` observations without requiring them to lie between the other grids. Do not require every individual translation, correction case or successive grid to improve monotonically: cancellation, interpolation phase and FP32 accumulation can make individual errors nonmonotonic. Save raw output for all grids before applying acceptance gates.

Passing establishes fixed-state classical electrostatic assembly, correction conventions and scoped mesh convergence for these small orthogonal cases. It does not establish arbitrary PME tolerance accuracy, triclinic MD, biomolecular force-field accuracy, energy conservation over long trajectories, NVT/NPT sampling or performance.

## Priority 2: constraints and Langevin statistics

### Finite-time unforced Langevin check

The frozen `ClassicalThermalConformanceTests` fixture uses 1,024 neutral noninteracting particles of mass 12 Da, eight declared seeds and 32 steps per seed, with observations at steps 1, 8 and 32. Temperature is 300 K, timestep is 1/256 ps, friction is 32/ps, and common initial velocity is `(0.25, -0.125, 0.0625)` nm/ps. Each horizon contains 8,192 particle/seed samples; different horizons are not pooled. Keep particles well separated and check that the cutoff/neighbor geometry does not introduce rejected steps. For this fixture the exact finite-time expectation is

`E[v_n] = exp(-gamma*n*dt) * v_0`

`Var[v_n] = (kB*T/m) * (1 - exp(-2*gamma*n*dt))`.

Mean and known-mean central second moments use Gaussian concentration bounds declared below, including cross-component moments and the full statistical comparison budget. This exercises the actual [native stochastic update](../../Sources/NumiVivoShaders/Resources/NumiVivoMDKernels.metal), line 29, and its time/seed indexing. It does not qualify a force-dependent configurational distribution.

### Constrained velocities and coupled constraints

At a fixed valid geometry, let `M` be the diagonal Cartesian mass matrix and let the rows of `G` be the independently constructed distance-constraint gradients. With independent constraints and no removal of center-of-mass momentum, the projected thermal Gaussian has covariance

`Cov(v) = kB*T * [M^-1 - M^-1*G^T*(G*M^-1*G^T)^-1*G*M^-1]`.

The test oracle constructs this small matrix directly, without calling the production constraint projector or CPU velocity initializer. The frozen fixture uses 128 rigid dimers of masses 4 and 16 Da, with relative coordinates `(0,0,0)` and `(0.125,0.09375,0)` nm, and 32 distinct `thermalize` seeds at fixed geometry. That gives 4,096 molecular samples. Each isolated dimer has five kinetic degrees of freedom: `2*K/(kB*T)` has the chi-square distribution with five degrees of freedom; center-of-mass kinetic energy has expectation `3*kB*T/2`, and rotational energy has expectation `kB*T`. Independently recompute distance and tangent-velocity residuals from readback.

The noncollinear rigid triatomic uses masses 12, 16 and 20 Da at `(0,0,0)`, `(0.125,0,0)` and `(0.03125,0.125,0.0625)` nm, with all three distance constraints, 128 copies and 32 seeds. Its six kinetic degrees of freedom exercise coupled iteration. The production default 32-iteration budget remains fixed; required relative distance residual is at most `1e-6` and tangent speed at most `2e-5` nm/ps. The reference geometry is checked against the initial accepted coordinates. A failure to converge must remain a failure of that supported fixture, not disappear through an undocumented iteration increase. Redundant or degenerate constraint rows need a distinct rank/compatibility decision and must not be assigned the independent-row chi-square expectation.

### Frozen statistical bounds

For the whole thermal family, declare `alpha = 1e-4`, comparison budget `B = 512`, and `t = log(2*B/alpha)`. The implemented counters require exactly 140 statistical comparisons (27 OU, 43 dimer and 70 triatomic) and reject a budget overrun. This union bound does not require independence between different comparisons. Deterministic state, oracle and constraint checks are separate.

For `N` independent Gaussian vector samples with known covariance `C`, the absolute mean allowance is `sqrt(2*Cii*t/N)` and the known-mean second-moment allowance is `sqrt(2*(Cii*Cjj+Cij^2)*t/N) + (sqrt(Cii*Cjj)+abs(Cij))*t/N`. Explicit FP32 allowances add `1e-5*sqrt(Cii)` and `1e-5*sqrt(Cii*Cjj)` respectively. Samples are centered about the analytical mean, not a fitted sample mean.

For total, molecular COM, rotational and whole-system COM kinetic energies, use the declared chi-square degrees of freedom. The mean lower/upper deviations are `2*sqrt(nu*t/N)` and `2*sqrt(nu*t/N)+2*t/N`, plus `1e-5*nu`. CDF values at `0.5*nu`, `nu` and `2*nu` use the Hoeffding allowance `sqrt(t/(2*N))+1e-5`. Whole-system COM has three degrees of freedom and 32 independent seed samples, catching unintended global COM removal. Sample kinetic variances are retained as observations without a separate normal-approximation acceptance gate.

Native `thermalize` samples and projects at fixed position ([runtime](../../Sources/NumiVivoKit/MD/VivoMDMetalRuntime.swift), lines 234–257). Passing its covariance test alone does not qualify constrained NVT trajectories. A subsequent fixed-length `dt` versus `dt/2` trajectory fixture should measure constraint residuals, center-of-mass and rotational energy statistics, and timestep bias with independent seeds and correlation-aware uncertainty. Long-run constrained sampling remains unqualified until that evidence exists.

## Priority 3: NPT proposal correctness before equilibrium claims

[pressureMove](../../Sources/NumiVivoKit/MD/VivoMDMetalRuntime.swift), lines 328–381, proposes a symmetric change in log volume and translates each molecular component's center while preserving its internal geometry. Its acceptance score at lines 366–368 contains `(componentCount + 1) * deltaLogVolume`. The extra one belongs to a proposal symmetric in log volume; it must not be replaced with an atom count or dropped as though the proposal were symmetric in volume.

A bounded first gate is 128 proposals for eight ideal monatomic particles and 128 for eight rigid dimers. Set every charge and Lennard-Jones coefficient to zero and use the cutoff profile. The frozen test declares temperature 300 K, pressure 150 bar, timestep 0.125 ps, initial side 2 nm, maximum log-volume step 0.25, and one fixed shared seed. Dimers have masses 1 and 3 Da and length 0.125 nm, including boundary-straddling molecules. Both systems contain eight molecular centers; their scalar volume and acceptance traces must agree despite their different atom counts.

Use zero initial velocities and zero Langevin friction to make the ideal dynamics stationary in exact arithmetic. The native runtime still executes FP32 periodic wrapping and constraint projection before each pressure proposal, so after a cell change the complete step may introduce rounding-scale coordinate or velocity changes. Save native pre-dynamics and post-step checkpoints, record exact-equality diagnostics, and compare accepted molecular-center geometry and rejected pre-step geometry with the predeclared physical tolerances. Save analytically proposed coordinates and cells for every proposal; the internal native pre-proposal/proposed buffers are not exposed by the public runtime. Accepted/rejected cell ownership and binary-step clock accounting remain exact. These are conditional transition checks, not equilibrium-volume or random-number-distribution qualification.

For each certificate, independently compute `V_proposed = V_before * exp(delta)` and

`score = -P * c * V_before * (exp(delta) - 1)/(kB*T) + (Nmol + 1)*delta`,

where `c = 0.0602214076` converts `bar * nm^3` to kJ/mol and `Nmol` is the number of independently translated molecular components. The oracle derives this conversion and `kB` independently from SI constants. Check the score within `1e-9`, accepted cell/volume within `1e-12` relative, position error within `1e-5` nm, COM error within `5e-6` nm, internal-length error within `6.25e-7` nm, internal-vector change within `2e-6` nm, and velocity magnitude/change within `5e-5` nm/ps. Potential energy remains within `1e-12` kJ/mol and kinetic energy within `1e-6` kJ/mol. Rejected cells remain exactly unchanged; full-step coordinate comparisons retain the rounding distinction above. A rejected certificate's `volumeAfterNM3` is the retained volume, not the proposed volume; derive the latter from `delta`. Require both accepted and rejected proposals and acceptance of every favorable proposal, keeping invalid-cell geometric rejection separate from ordinary Metropolis rejection.

This establishes a transition-score and restoration contract. The ideal-gas equilibrium volume density is a stronger expectation:

`p(V) proportional to V^Nmol * exp[-P*c*V/(kB*T)]`.

Without a lower-volume bound this is a Gamma distribution with shape `Nmol + 1`, scale `kB*T/(P*c)`, mean `(Nmol + 1)*kB*T/(P*c)` and variance `(Nmol + 1)*(kB*T/(P*c))^2`. The runtime's cutoff and neighbor-radius cell contract imposes a lower-volume boundary. For a cubic cell with neighbor lists, the side must be strictly greater than `2*(cutoff + skin)` ([radius validation](../../Sources/NumiVivoKit/MD/VivoMDMetalABI.swift), lines 137–154). A distribution test must compare with the correctly truncated density, or explicitly bound the excluded probability before using the untruncated moments.

A short proposal test is not an NPT equilibrium qualification. That requires a separately frozen burn-in, trajectory length, independent-seed set and autocorrelation-aware effective sample size, with enough support for the stated confidence interval. Positive-friction dynamics and interacting constrained systems are further scopes. PME-at-two-volumes agreement helps validate electrostatic volume response, but does not by itself establish detailed balance or equilibrium sampling of the coupled PME/NPT system.

## Evidence required when executing this plan

Record the exact candidate commit and source archive, native host/toolchain, selected tests, seed set, configuration, fixed native-work bounds and unmodified predeclared thresholds. Save reference convergence errors, per-case energy and force errors, pooled mesh trend, checkpoint-preservation checks, stochastic sample counts/uncertainty and any rejection certificates. Missing native hardware or an unconverged reference is not a pass. Keep build success, successful execution and scientific acceptance separate, and publish measured results in a later evidence record rather than rewriting this source audit as though it had measured them.
