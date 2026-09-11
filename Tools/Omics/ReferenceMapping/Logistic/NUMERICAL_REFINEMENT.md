# Numerical convergence refinement, 2026-09-11

The first complete native run used gradient tolerance 1e-7. All four fits/maps
and native reconstructions completed. The independent checker stopped on human1
because maximum probability difference from its tighter scikit-learn fit was
0.00014894968297529676, exceeding the previously fixed 0.0001 acceptance bound.
That assertion precedes reading query outcome labels or calculating new accuracy
metrics. All four first-run native archives and the failed check are retained.

Refit every original fold at gradient tolerance 1e-8 with the same executable,
objective, penalty, class weights, standardization, count sources, labels and
query plans. This tightens numerical convergence; the probability comparison
bound remains 0.0001. No donor/class/feature selection, biological parameter or
outcome criterion changes. The original protocol and first input freeze remain
intact; a new input freeze binds the single tolerance change and reused gzip
sources before fitting. The product option's default remains 1e-7 and its
reported gradient criterion is explicit. The qualified comparison uses 1e-8.
