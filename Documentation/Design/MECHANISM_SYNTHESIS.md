# Abstract mechanism synthesis

## Purpose

NumiVivo converts a behavioral program into candidate abstract molecular mechanisms. This stage selects compatible functional parts and predicts computational performance envelopes. It does not produce nucleotide sequences, laboratory procedures, delivery instructions, or clinical recommendations.

The synthesizer operates on two documents:

- a **problem** containing required functional slots and global constraints;
- a **library** containing sequence-free mechanism records and evidence envelopes.

## Functional roles

A part has one role:

- sensor;
- transducer;
- logic;
- temporal filter;
- memory;
- effector;
- communication;
- containment;
- shutdown;
- monitor.

A complete design assigns at most one part to each slot. Optional slots may remain empty. A future multi-part slot format must be versioned rather than overloading the current one-part semantics.

## Part contract

A part declares:

- stable identifier;
- role;
- accepted abstract signals;
- produced abstract signals;
- supported hosts;
- supported delivery modes;
- dependencies;
- pairwise incompatibilities;
- orthogonality group;
- evidence tier;
- reversibility class;
- independent-control status;
- context-insulation status;
- resource-buffering status;
- performance envelope.

Empty host or delivery sets mean that the library author has not imposed a restriction. They do not establish broad biological compatibility.

## Performance envelope

Each part contains bounded model estimates for:

- payload size;
- cellular burden;
- latency;
- leakage probability;
- dynamic range;
- specificity;
- robustness;
- relative uncertainty.

Candidate aggregation uses explicit rules:

```text
payload       = sum(part payload)
burden        = sum(part burden)
latency       = sum(part latency)
leakage       = 1 - product(1 - part leakage)
dynamic range = minimum(part dynamic range)
specificity   = minimum(part specificity)
robustness    = minimum(part robustness)
uncertainty   = root-sum-square(part uncertainty)
```

These rules are computational assumptions. A MechanismPack must record the aggregation version so later models can use different composition laws without changing historical results.

## Evidence

Evidence tiers are ordered from strongest to weakest:

1. observed in the declared target context;
2. observed in a related context;
3. calibrated;
4. inferred;
5. assumed;
6. hypothetical.

The weakest selected part becomes the candidate's evidence tier. Evidence does not become stronger when several weak records are combined.

## Reversibility

Reversibility classes are:

- reversible;
- resettable;
- conditionally irreversible;
- irreversible.

The least reversible selected part determines the candidate class. Slots may impose a maximum allowed class. An independently controlled shutdown slot must be reversible or resettable in the current synthesis contract.

## Compatibility

Before search, the compiler filters parts by:

- role;
- host;
- delivery mode;
- required input;
- required output;
- evidence limit;
- reversibility limit;
- independent-control requirement;
- context-insulation requirement;
- resource-buffering requirement;
- required tags.

During search, it rejects:

- explicit pairwise incompatibility;
- duplicate orthogonality groups when required;
- unsatisfied dependencies;
- missing independent shutdown;
- missing monitor;
- global performance limits.

Dependencies must refer to another selected part identifier. A dependency does not automatically add a part to the candidate.

## Search order

Slots are ordered by increasing candidate count. This exposes contradictions early and reduces the branch count.

For each partial assignment, the synthesizer computes aggregate metrics and rejects the branch as soon as it exceeds any monotonic upper bound. It also prunes a branch when its current objective cannot improve on the worst retained candidate after the result set is full.

Search is bounded by:

- maximum visited nodes;
- maximum returned solutions;
- document parser limits;
- part and slot collection limits.

Exhausting the search budget returns valid candidates found so far and an explicit warning. It does not present them as a proven global optimum.

## Pareto filtering

A candidate dominates another candidate only when it is no worse on all tracked dimensions and strictly better on at least one. Dominated candidates are removed before the maximum solution count is applied.

The remaining candidates are ordered by a configurable scalar objective. Objective weights affect ordering but do not replace hard constraints or Pareto filtering.

## Candidate identity

Candidate identity is a SHA-256 digest over:

- synthesis format version;
- problem identifier;
- host;
- delivery mode;
- sorted slot-to-part assignments.

The identity does not include an optimizer iteration count or search order. Identical functional assignments therefore retain the same identity across deterministic search refactoring, provided format semantics remain unchanged.

## Diagnostics

The result contains rejection counts for:

- host incompatibility;
- delivery incompatibility;
- signal incompatibility;
- weak evidence;
- reversibility mismatch;
- missing independent control;
- incompatible parts;
- duplicate orthogonality group;
- missing dependency;
- payload, burden, latency, leakage, uncertainty, specificity, robustness, or dynamic-range violation;
- missing safety role;
- search-budget exhaustion.

Rejection counts support model debugging but are not probabilities.

## C ABI

The public ABI accepts bounded UTF-8 JSON problem and library documents. It returns:

- candidate JSON;
- diagnostic JSON;
- a typed status code.

Returned buffers use the same caller-release contract as the main compiler ABI. The ABI does not execute scripts, access the network, load dynamic libraries, or resolve external sequence databases.

## Integration with VivoProgram

The intended compiler sequence is:

```text
VivoProgram
  -> typed behavior IR
  -> functional mechanism slots
  -> constrained mechanism synthesis
  -> selected abstract mechanism graph
  -> reaction and regulation IR
  -> ProgramPack + MechanismPack
```

A selected candidate remains abstract until an independent downstream system resolves it into a physical design under appropriate controls. NumiVivo records that boundary explicitly.

## Safety boundary

A mechanism candidate is a model result. It does not establish:

- physical realizability;
- host compatibility;
- absence of off-target behavior;
- delivery feasibility;
- containment effectiveness;
- therapeutic efficacy;
- safety in an organism;
- regulatory acceptability.

No simulation or synthesis score authorizes physical implementation.
