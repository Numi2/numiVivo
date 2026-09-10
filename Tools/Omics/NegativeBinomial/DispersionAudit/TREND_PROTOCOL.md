# Trend-stage decomposition

Declared after the first stage audit results exposed different curves despite
equal Hagai prior variance, before fitting any of the following diagnostic
curves. Apply the same set to all twenty splits; preserve every warning/failure.

Use pinned DESeq2's `parametricDispersionFit` on (1) its original gene-wise
estimates at its original admission threshold 1e-6, (2) native gene-wise
estimates on the exact native trend reference set, (3) DESeq2 estimates on
that same native set, and (4) native estimates on the native reference set
restricted to at least 1e-6. The original native/R estimates, curves, tested
families and calls stay unchanged. This is an input/algorithm diagnostic,
not a candidate inference policy or a new false-discovery experiment.

Record reference-set sizes, coefficients and maximum curve-relative differences
over the observed means. Case (1) must reproduce the stored DESeq2 curve. Do not
assume case (2) matches native: the two implementations' iteration and admission
rules remain independently authored and any discrepancy is part of the result.
