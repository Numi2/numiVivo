# Explicit support-work sensitivity

The complete first pass exhausted the frozen one-million-point per-observation
budget in twenty gene/arm fits: nine Hagai treatment genes under each trend,
and one Crowell excitatory-neuron gene under each trend. These remain failed
primary stages. No mean, dispersion, QL scale, prior, tolerance or feature family
is changed to remove them.

Before evaluating their missing moments, rerun all twenty exact failed gene/arm
inputs through the same native adjustedResiduals owner, with its public
maximumTerms explicitly set to 10,000,000. Keep relativeTolerance=1e-10 and
all original counts, means, design and supplied scales. Record each result,
work count and failure. This is a computational work-limit sensitivity, not a
new method selected by statistical outcomes, and does not overwrite primary
failures or change the default budget. Independently check every resulting
moment with PMF summation, retaining any numeric disagreements.
