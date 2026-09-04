# Exact and hybrid stochastic execution

## Objective

NumiVivo represents molecular systems whose relevant populations range from individual events to effectively continuous concentrations. A single numerical method is not authoritative across that entire range.

The stochastic stack therefore contains:

- exact direct-method stochastic simulation for low-count and critical components;
- bounded tau-leaping for event-dense discrete components;
- deterministic RK2 for high-count continuous components;
- spatial split execution for components coupled to transport fields;
- a graph planner that prevents overlapping numerical authority.

## Reaction-component ownership

Before selecting a method, the planner builds a reaction hypergraph:

- every reaction is a node;
- every dynamic species connects all reactions that read or write it;
- connected reactions form one authority cohort.

This prevents an exact reaction and a leap reaction from consuming the same species independently during one transaction.

A component can be split only through an explicit operator-splitting contract that introduces a publication boundary and an error measure. The default planner does not split connected molecular pools.

## Exact direct method

For state `x`, reactions `j`, and propensities `a_j(x)`, the direct method evaluates:

```text
a_0 = sum_j a_j
τ = -ln(u_1) / a_0
choose μ such that cumulative propensity crosses u_2 a_0
x <- x + ν_μ
```

The Metal implementation assigns one independent lane to one GPU thread. Each lane advances through a bounded event loop until it reaches the target time or event budget.

The exact runtime uses UInt32 authoritative counts and signed integer stoichiometric changes. Before applying an event it verifies:

- propensity is finite and non-negative;
- reactant counts are available;
- signed updates cannot underflow or overflow;
- event time is finite and monotonic;
- the per-step event limit is not exceeded.

Exceeding the event limit does not truncate the trajectory. The candidate is rejected or retried with a smaller time window.

## Counter-based random streams

Random identity is derived from immutable identifiers rather than scheduling order. The conceptual namespace is:

```text
program fingerprint
experiment fingerprint
replicate seed
lane identity
reaction component identity
logical step
substep or retry identity
event counter
sample lane
```

A checkpoint stores the random-stream version and all counters required to continue the same stream. A migration may preserve a stream only when the owning component and stochastic semantics remain compatible.

## Propensity classes

The exact model supports explicit propensity classes with validated parameters and reactant references. The current abstract model contains:

- zero-order source;
- first-order reaction;
- second-order reaction with distinct reactants;
- second-order reaction using two copies of one reactant;
- Hill activation;
- Hill repression.

A rate-law implementation is accepted only when its required reactants and parameters are present. Unsupported laws produce a compile diagnostic rather than falling back to mass action.

## Tau-leaping

Tau-leaping estimates the number of firings of each reaction over a bounded interval. NumiVivo adds two protections:

1. Critical reactions use a bounded event treatment rather than an unconstrained Poisson draw.
2. Competing firings are scaled against shared reactant availability before state updates are published.

A leap candidate is rejected when it produces invalid counts, non-finite propensities, excessive event volume, or a monitor violation.

The leap error target and maximum expected firings are part of the hybrid policy and plan certificate.

## Deterministic execution

High-count components may use continuous FP32 state when:

- population statistics are above the deterministic-entry threshold;
- expected event volume exceeds the stochastic efficiency envelope;
- the component contains no reaction requiring exact event semantics;
- a declared continuous approximation is permitted by the requested fidelity.

Deterministic entry and exit thresholds differ to prevent method oscillation around one count boundary.

## Spatial execution

A reaction component with spatial transport dependence is assigned to the spatial-split authority. Reaction and transport stages are ordered explicitly. A method change cannot silently discard voxel topology, boundary conditions, velocity fields, or diffusion coefficients.

## Planner inputs

Each reaction descriptor supplies:

- reaction index;
- dynamic species indices;
- critical flag;
- spatial flag;
- minimum and median reactant counts;
- expected firings per proposed step;
- maximum propensity;
- stiffness ratio.

The planner uses these summaries to choose an authority and a conservative maximum step. The plan records rationale for every cohort.

## Hysteresis

Method changes have separate entry and exit thresholds:

```text
exact entry count < exact exit count
continuous exit count < continuous entry count
```

An exact cohort remains exact until it clears the wider exit threshold. A continuous cohort remains continuous until it falls below the lower exit threshold. This reduces repeated representation conversion.

## State migration

### Continuous to discrete

The migration kernel requires every value to be:

- finite;
- non-negative;
- no larger than the UInt32 count limit;
- within the declared tolerance of an integer.

A materially fractional population is not rounded silently. The migration is rejected.

### Discrete to continuous

UInt32 values above the exact binary32 integer range are reported. Policy may keep the cohort discrete, select a wider continuous representation in a later format, or reject migration. Silent loss of integer identity is prohibited.

### Preserved state

Random counters and delayed queues are preserved only when the migration plan says their semantics remain stable. Spatial and tissue topology is rebuilt when entering a representation that requires it.

## Hybrid plan integrity

A hybrid plan contains:

- species count;
- independent cohorts;
- authority mode per cohort;
- reaction and species membership;
- recommended maximum step;
- leap error target or exact-event budget;
- previous mode;
- migration requirement;
- human-readable rationale;
- deterministic fingerprint.

Every species owned by a reaction component appears in exactly one cohort. Species not participating in the reaction graph are listed separately.

## Failure modes

The runtime distinguishes:

- invalid propensity;
- count underflow;
- count overflow;
- non-finite event time;
- exact-event budget exhaustion;
- leap-event capacity exhaustion;
- representation-conversion rejection;
- random-stream incompatibility;
- delayed-queue incompatibility;
- spatial-topology migration failure.

Failures retain the previous committed state.

## Interpretation

Exact stochastic execution is exact with respect to the declared Markov reaction model and its propensity implementation. It does not establish that the declared reactions, rates, compartments, or host context are biologically complete or experimentally correct.
