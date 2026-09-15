# Packed Metal negative-binomial objective experiment

## Scope

`VivoMetalNegativeBinomialLikelihood.evaluateBatch(_:)` is an internal,
unpromoted arithmetic primitive. It packs independent, fixed-dispersion,
FP32-domain likelihood objectives into shared Metal dispatches, while keeping
each request's input order and compensated FP64 CPU reduction separate.

This is deliberately narrower than a cohort backend. It neither changes the
default CPU likelihood nor submits a biological or statistical result from the
GPU. The existing opt-in `metalFP32` fitter remains experimental. Its v2
execution profile now also records whether the current FP32 Metal line-search
decision agrees with the exact FP64 likelihood for each candidate in a cohort
run. That audit does not alter the chosen candidate, recover failed fits, or
describe this packed primitive.

## Contract

Each `BatchRequest` has its own count vector, mean vector, and fixed
dispersion. Counts, dimensions, means, and `dispersion * mean` must meet the
same bounded FP32 contract as the scalar objective. A packed dispatch writes
one FP32 term per observation; the host reduces the terms in request and
observation order with a separate compensation term for each request.

The packed kernel has a per-observation dispersion buffer because independent
genes can have different fixed dispersions. It is therefore not valid to add
their terms into one likelihood or to reuse a result for another gene.

## Deliberate non-claims

This primitive is not yet connected to the serial coefficient/dispersion
solver. Consequently it is not evidence of end-to-end speed, numerical
coverage, line-search behavior, fit recovery, statistical agreement, or
biological validity. In particular, it does not alter the Kang benchmark
receipt or make the Metal option promotable.

Before a cohort scheduler can use it, it must preserve deterministic gene
ordering, submit only independent fixed-dispersion requests, retain the CPU
acceptance owner, and be remeasured on the declared supported domain against
the CPU-successful fits. The focused Metal test compares packed values with
both the scalar Metal objective and an FP64 CPU oracle; it is an arithmetic
guard, not a performance qualification.

The v2 profile audit is deliberately a failure-isolation step before changing
the solver: it records both directions of decision disagreement and the maximum
normalized margin difference under the current objective tolerances. It is not
evidence that FP32 line-search decisions explain every current failure, and it
does not replace a profiler or a complete same-domain CPU coverage/performance
comparison.
