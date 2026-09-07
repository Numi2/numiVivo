# Native constant-pH NCMC

`VivoConstantPHNCMCMetalFactory` connects the existing constant-pH sampler to the
existing Metal MD evaluator and integrator. No second classical, QM/MM,
polarization, or virtual-site force implementation is introduced.

```swift
let setup = try VivoConstantPHNCMCMetalFactory.make(
    specifications: preparedStateSpecifications,
    samplingTemperatureK: 300,
    configuration: .init(
        schedule: .linear(intervals: 16, propagationStepsPerLambda: 4),
        timeStepPS: 0.0005
    )
)
let sampler = try VivoConstantPHNCMC(
    configuration: .init(base: semigrandConfiguration),
    initialPhysicalState: acceptedPhysicalState,
    executableStates: setup.executableStates,
    switchEngine: setup.switchEngine
)
let result = try await sampler.run()
```

The caller supplies prepared, parameterized endpoint Hamiltonians with the same
physical particle identities, masses, constraints, dependent-site constructions,
and fixed cell. Proton-changing preparation must use the explicit latent-proton
union topology and its reservoir calibration; the factory does not create pKa
parameters. State identifiers and complete endpoint execution fingerprints are
bound to the switch engine and checked before every proposal.

## Hamiltonian and proposal

The alchemical path is `U(l) = (1-l) U_A + l U_B`, where both endpoints include
all retained classical terms and any electronic/polarization provider. Each
endpoint is evaluated through `VivoMDMetalRuntime.evaluateHamiltonian(at:)`.
Endpoint forces are already redistributed to physical particles; the invariant
mass/constraint/site carrier does not add their potential a second time.

Every lambda interval is **perturb to midpoint, propagate at midpoint, perturb
to endpoint**. The reverse edge uses the mirrored nonuniform grid
`1 - reverse(lambdas)`. Thus reversing the endpoint states and momenta reverses
the actual operation order, not just the labels. A plain repeated
perturb-then-propagate schedule does not have that property.

Fixed-lambda propagation is NVE. Its measured shadow work is `delta(U+K)`, not
zero merely because velocity Verlet is symplectic. The total protocol plus
shadow work is checked against the complete endpoint Hamiltonian difference.
All work is kJ/mol; `Da*(nm/ps)^2` is kJ/mol. The outer Metropolis-Hastings rule
adds the explicit semigrand reservoir difference and neighbor-proposal ratio.
On rejection, pre-switch positions are restored and pre-switch velocities are
reversed. On acceptance, the proposed phase-space state is retained.

Version 2 rejects raw stochastic propagation path ratios: a thermostat needs an
explicit heat/action convention to prevent double counting. No such quantity is
inferred from endpoint energies. Nonfinite work, changed cells/manifolds,
unconverged endpoint providers, cancellation, or capacity exhaustion fail before
publishing the candidate. A completed progress callback observes a committed
checkpoint; a callback error does not undo that already published transaction.
Checkpoint validation reconstructs the discrete state/RNG/acceptance history,
but does not replace re-evaluation of physical trajectory evidence.

## Execution and scope

Endpoint runtimes are reused within a proposal. Only exact-geometry endpoint
values may be cached across lambda changes. At most two endpoint runtimes and
one midpoint runtime are live. Work budgets count endpoint evaluations, not
only accepted steps. Independent positions and velocities must explicitly
round-trip through FP32; no unaccounted coordinate rounding is accepted.

The target uses the declared invariant masses. Metal trajectory arithmetic and
constraint tolerance remain finite-precision approximations: energy accounting
alone does not establish numerical reversibility, equilibrium convergence, or
chemical accuracy. NVE endpoint propagation is useful for tests, but production
configurational sampling must also establish the intended canonical ensemble.
Linear endpoint mixing is **not** soft-core alchemy; singular or unconverged
endpoint evaluations are errors. Massive Drude variables, changing constraints,
cell moves, adaptive mass manifolds, and stochastic NCMC are not supported here.

The NCMC configuration/checkpoint/result and lambda schedule/result schemas are
v2. Old v1 evidence is not silently reinterpreted after changing proposal order
and momentum semantics. `VivoNCMCProtocolKernel.switchEngine` is now throwing and
binds its schedule/work limits into its returned fingerprint.

`NCMCProtocolKernelTests` checks nonuniform reverse schedules and antisymmetric
protocol/shadow work. `ConstantPHNCMCTests` checks seed-zero replay, rejection
momenta, work limits and tampering. `ConstantPHNCMCMetalTests` exercises actual
Metal propagation, complete external-provider participation, energy telescoping,
round trips, checkpoint continuation, and admission/rounding rejection. These are
software fixtures, not a protein titration benchmark.

Method references: Nilmeier et al., PNAS 2011, DOI 10.1073/pnas.1106094108;
Gill et al., JPCB 2018, DOI 10.1021/acs.jpcb.7b11820. The implemented profile is
deterministic NVE NCMC, not the stochastic protocol of every cited application.
