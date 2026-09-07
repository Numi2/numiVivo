# Nuclear statistics, global uncertainty and reactive delta sampling

These are reusable extensions of the existing molecular Hamiltonian and workflow
interfaces. None creates another QM/MM electronic or classical force field.
Each extension has a distinct scientific scope: one-dimensional barrier models,
equilibrium nuclear statistics, nonlinear propagation of supplied joint evidence,
and accelerated equilibrium sampling. They do not certify protein chemistry.

## 1. Quantum-nuclear calculations

`VivoBarrierTunnelling.calculate` evaluates Wigner or asymmetric Eckart corrections.
The request binds temperature, a positive imaginary-frequency magnitude in cm^-1,
positive forward/reverse barriers in kJ/mol, one stated energy/ZPE convention,
the Hamiltonian fingerprint and stationary-point evidence. Wigner is admitted
only for h*nu/(kBT) <= 1. Eckart uses analytic channel transmission, logarithmic
thermal integration, two-resolution convergence, and a bounded high-energy tail.
Closed channels have zero transmission. A shallow barrier uses the cosine
analytic continuation, not the absolute value of a negative channel argument.
Forward/reverse barrier exchange preserves the correction factor. Results retain
model warnings, convergence statistics and the logarithmic factor.

`VivoNuclearRateComposition` is a separate explicit model composition. Temperature
and Hamiltonian must match; an already quantum-corrected rate is rejected. The
caller must state the separability assumption. It does not overwrite the
classical reactive-flux transmission coefficient or silently amend a kinetic pack.
Numerical convergence of an Eckart integral does not establish the validity of a
one-dimensional tunnelling approximation. In particular, these are not a
multidimensional instanton or a quantum dynamical recrossing calculation.

`VivoRingPolymerSampling` samples the distinguishable-nucleus, finite-bead primitive
path integral at fixed temperature and cell. It uses the same complete potential
at every bead, an exact normal-mode free-spring propagation, physical-potential
kicks, refreshed momenta and a full-Hamiltonian Metropolis decision. With P beads,

```
omegaP = P*R*T/hbar_molar
UP = sum_b [U(q_b) + sum_i m_i*omegaP^2*|q_bi-q_(b+1)i|^2/2]
probability(q) proportional to exp[-UP/(P*R*T)]
```

Energies are kJ/mol, coordinates nm, masses Da, sampler velocities nm/ps. The
normal modes avoid numerically integrating the stiff free-ring springs. Rejected
states remain in the chain. The primitive energy estimator, centroid, action,
acceptance, force-call accounting and exact restart state are retained. The
reduced configurational potential includes the `-3P/2 * sum(log(mass))` mass
normalization for same-P, same-T isotope reweighting. It is not a standalone
absolute partition function. Bead-number convergence, equilibration, correlation
and cross-isotope overlap require their own calculations. Existing multistate
MBAR accepts the resulting cross-evaluated reduced potentials.

The executable force contract is `VivoNuclearPotential`. It contains an explicit
atom-to-particle layout, masses, Hamiltonian identity and coordinate precision.
`VivoNuclearMetalSpecification` constructs classical, fixed QM/MM or adaptive
QM/MM potentials through the existing force factory and `VivoMDMetalRuntime`.
Dependent sites are reconstructed by that authority; returned forces are physical
forces with dependent-site contributions already redistributed.

The current native runtime stores coordinates in FP32. Its adapter therefore
explicitly declares the numerical potential `U_FP32(round(q))`; it does not claim
a continuous FP64 force authority or weaken the native checkpoint validator.
Metropolis acceptance uses that same energy definition. Native projection and
continuous FP64 reference profiles are not interchangeable checkpoint identities.
For a periodic system, all beads and physical atoms of a molecule receive the
same lattice translation when choosing its coordinate chart. No independent atom
wrapping is used. This is a zero-winding molecular path model, not exchange or
winding-number sampling. Constraints, massive Drude particles and cell moves are
not accepted by this independent-coordinate sampler.

## 2. Global joint uncertainty

`VivoGlobalJointUncertainty.propagate` evaluates complete joint parameter draws,
not independent samples of marginal error bars. Every draw carries a model ID,
log weight, declared dependence block and a complete vector in the declared
parameter order. Model probabilities are normalized independently of the number
of draws supplied for each model. All model probabilities and their provenance
are explicit inputs. An optional joint Gaussian constructor accepts a supplied
full positive-semidefinite covariance, including singular covariances; it never
adds jitter that creates unsupported independent uncertainty.

Each full draw is propagated through the nonlinear forward calculation. Outputs
include the full weighted ensemble, mean, population variance, weighted quantiles,
within-model and between-model variance, and cross-observable covariance. These
are distribution summaries, not Monte Carlo standard errors of the mean. An
undefined observable retains its probability of being defined and conditional
summary; it is not silently zero-filled or deleted from unrelated observables.
Numerical forward failures abort the result instead of renormalizing over only
successful draws. Weight and declared-block effective sample sizes do not claim
measured autocorrelation, independent trajectories or confidence-interval coverage.

`VivoGlobalKineticUncertainty` provides the typed chemistry implementation. It
first validates the existing state-specific replicated rates, exchanges and
chemical context, then applies joint changes to relative state free energies,
initial-population logits, pathway rates, electronic barrier shifts, classical
transmission factors, explicit nuclear factors and state-exchange rates. It
recalculates populations and evaluates the resulting generator with the existing
uniformization solver. Survival, reacted probability, hazards and state-resolved
probabilities retain their correlations at every observation time. Classical
transmission remains in its physical interval; invalid draws are rejected, not
clipped. Underflowed hazards remain undefined.

This supplies global nonlinear propagation, not automatic Bayesian inference of
missing experimental evidence. Baseline assumed populations, exchange rates,
transmission, state free energies and nuclear separability remain visible in the
assumptions. Unresolved components remain unresolved. An unspecified chemical
state or alternative pathway is not made present by adding a covariance matrix.
The result does not automatically pass the kinetic-pack qualification gate.

## 3. Reactive delta surrogate and accelerated equilibrium sampling

`VivoReactiveDeltaSurrogate.train` learns the difference between a supplied baseline
and authoritative complete potential. It fits energies and analytic forces jointly
using bounded Gaussian radial-basis features of explicitly mapped pair distances.
There are no fixed bond predicates in the features, so the declared local domain
can include changing bond distances. The model is specific to its mapped atoms,
Hamiltonians and training domain, not a transferable chemical force field.

Train and held-out source groups must be disjoint. Basis centers come only from
training data; duplicate geometry leakage is rejected. Bootstrap committee members
resample complete training groups. Held-out energy/force tolerances qualify the
model; distance-to-training support and committee disagreement determine whether
it may propose a trajectory. Committee disagreement is a diagnostic, not a
calibrated posterior uncertainty. Acquisitions must be grouped by source trajectory
or chain rather than by individual correlated frame.

`VivoReactiveEnergyForceBackend.metal` is an actual batched Metal implementation of
the radial-basis energy and feature derivatives. It uses FP32 inference and mapped
force reconstruction, without float atomics. FP64 fitting and authoritative energy
accounting remain separate. CPU and GPU inference profiles are fingerprinted and
are not silently exchanged during restart. Numerical acceleration is implemented;
a throughput or speedup claim requires a measured workload and reference timing.

`VivoReactiveSurrogateSampling` uses a frozen qualified model to generate reversible
baseline-plus-delta proposals and evaluates the authoritative complete energy at
proposal endpoints. Acceptance includes the full authoritative energy and kinetic
change. It therefore does not treat the learned force as chemical authority.
Every rejection remains an equilibrium sample. Domain failures create explicit
acquisition requests, not an unrecorded switch in force law halfway through a
trajectory.

To avoid restricting the accessible state space to the learned domain, a fixed,
state-independent positive probability selects a Gaussian random-walk proposal
accepted using the exact energy. Those moves can leave the surrogate support.
The mixture probability, step size and all RNG draws are checkpointed. The
surrogate-only path is not allowed to become the sole support-limiting kernel.
Whole-molecule periodic chart translations preserve proposal symmetry. There is
no selection of a kernel based on the current committee disagreement.

Authoritative endpoint labels, including rejected endpoints, are retained. A model
is frozen for an entire checkpointed epoch. Training a new model creates a new
model identity and requires a new epoch; it cannot rewrite past acceptance ratios.
This is equilibrium Monte Carlo sampling, not a physical MD clock. Its trajectories
must never be supplied as unbiased NVE shooting or used to calculate classical
recrossing, transport or residence times.

## Workflow and CLI use

The normal workflow registry exposes the following operations. Deterministic
operations reconstruct their numerical output. Native sampling operations validate
source identities, exact checkpoint transitions, stored accounting and proposal
RNG replay; cached validation is not a new GPU trajectory or chemical experiment.

| Operation | Output |
|---|---|
| `vivo.platform.barrier-tunnelling` | `correction` |
| `vivo.platform.global-kinetic-uncertainty` | `prediction` |
| `vivo.platform.reactive-surrogate-label` | `labels` |
| `vivo.platform.reactive-surrogate-train` | `model` |
| `vivo.platform.ring-polymer-sample` | `run`, `checkpoint` |
| `vivo.platform.reactive-surrogate-sample` | `run`, `checkpoint` |
| `vivo.platform.ring-polymer-convergence` | `assessment` |
| `vivo.platform.reactive-surrogate-coverage` | `assessment` |

Use `VivoRingPolymerWorkflowRequest`, `VivoReactiveSamplingWorkflowRequest` and
`VivoReactiveLabelWorkflowRequest` for native request construction. Potential
specifications are serializable configuration; executable closures are constructed
locally by the force factory, never deserialized. Resource admission accounts for
simultaneously resident classical systems, electronic scratch, meshes, labels,
beads, model inference and the requested number of force calls.

A runnable one-dimensional numerical fixture is provided:

```sh
swift build -c release
.build/release/numivivo workflow-plan Examples/workflows/barrier-tunnelling.json
.build/release/numivivo workflow-run Examples/workflows/barrier-tunnelling.json \
  --store /tmp/numivivo-nuclear-fixture --output /tmp/numivivo-nuclear-result.json
swift test -c release --filter 'BarrierTunnellingTests|RingPolymerSamplingTests|GlobalJointUncertaintyTests|GlobalKineticUncertaintyTests|ReactiveSurrogateTests|PrecisionMetalSamplingTests|PrecisionPlatformTests'
```

The Apple gate includes real native potential evaluation, ring-polymer restart,
Metal batched surrogate inference, label generation, training and authoritative
surrogate-Monte-Carlo execution. Synthetic fixtures test algorithmic contracts,
not protein chemistry. `Tools/PrecisionSampling/reference_checks.py` independently
checks the Eckart thermal integral with SciPy adaptive quadrature; it is a reference
tool, not a native-runtime dependency.

## Qualification campaigns

`VivoCorrelatedSamplingAnalysis` computes classical split-R-hat, a conservative
initial-positive-sequence autocorrelation effective sample size and Monte Carlo
standard error from equal-length independent chains. Constant but displaced
chains explicitly fail mixing. The autocorrelation work is bounded; a sequence
that remains positive at that bound is marked truncated and cannot pass. These
are bulk scalar diagnostics with at most ten million retained scalar values; they
are not rank-normalized tail diagnostics and do not prove that the discarded
interval removed initialization bias.

`VivoRingPolymerConvergence.assess` binds complete retained runs to one potential,
temperature, declared observable and predeclared selection protocol. Every bead
level requires distinct checkpoint starts and at least two chain identities. A
level must pass the chain diagnostics, and the requested number of consecutive
bead refinements must have an upper difference bound below the declared absolute
plus relative tolerance. The bound includes both levels' Monte Carlo standard
errors. Available observables are the primitive total-energy estimator, mean
potential energy, mass-weighted ring radius of gyration and a mapped centroid
distance. A passing synthetic harmonic campaign does not qualify a prepared
molecular Hamiltonian, isotope reweighting overlap or a quantum rate.

`VivoReactiveSurrogateCoverage.assess` evaluates a frozen model against grouped
authoritative labels from a separately fingerprinted campaign. Training-data
identity, duplicate geometries and reuse of held-out selection groups are rejected.
Each source group must independently meet descriptor/committee eligibility and
the frozen held-out energy/force limits, so one long correlated trajectory cannot
hide an unsupported region. Analysis is capped at 1,024 groups and the configured
primitive-work budget. The collection protocol and independence declaration remain
in the result because fingerprints cannot authenticate how evidence was generated.
Coverage still does not make learned forces authoritative or establish mixing,
throughput or speedup.

The next scientific runs must therefore supply multiple independently initialized
prepared-system chains at increasing bead counts and new grouped surrogate probes.
Until those evidence objects pass, the implementation remains executable but the
corresponding prepared system remains scientifically unqualified.

## Method references

- RMG-Py official Wigner and Eckart model documentation:
  https://reactionmechanismgenerator.github.io/RMG-Py/reference/kinetics/wigner.html
  and https://reactionmechanismgenerator.github.io/RMG-Py/reference/kinetics/eckart.html.
  NumiVivo explicitly normalizes the thermal integral in dimensionless energy and
  enforces the physical open-channel threshold.
- i-PI project, nuclear path-integral methods and implementation references:
  https://ipi-code.org/ . The finite-bead action above is independently checked
  against its analytic harmonic normal-mode distribution.
- Y. Nagai et al., *Self-learning Hybrid Monte Carlo: A First-principles Approach*,
  Physical Review B 102, 041124 (2020), https://arxiv.org/abs/1909.02255.
  This motivates authoritative endpoint acceptance rather than uncorrected learned
  forces; this implementation freezes its model during each sampling epoch.
