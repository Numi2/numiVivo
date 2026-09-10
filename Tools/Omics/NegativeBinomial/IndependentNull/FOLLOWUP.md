# Diagnostic availability follow-up

Declared after inspecting the original cohort 01 outputs, while the remaining
frozen native baseline was still running on 2026-09-10. This is a discovered
implementation repair, not a predeclared calibration analysis or a replacement
for `PROTOCOL.md`.

The baseline Wald/LRT owner unconditionally requires available Cook's distances.
An original source batch with a single donor produces a unit-leverage observation,
so those distances are undefined even when the treatment contrast is identified.
In cohort 01 this withholds 13,854 fitted genes, despite the request explicitly
having no influence threshold. Another 844 genes happen to pass the same guard,
making test availability depend on floating-point leverage rounding.

Require Cook's distances only when the caller requests an influence threshold.
Continue to retain nil diagnostics and withhold inference if such a requested
threshold cannot be evaluated. Preserve count, fit, dispersion, batch, donor,
normalization, support, boundary and BH rules. Do not invent influence values or
change the predeclared experimental-model defaults.

Add a focused singleton-batch regression covering both Wald and LRT, including
the explicit influence-policy failure. Retain the original 27 outputs. Run the
repaired Wald and LRT methods on all nine original cohorts as `native-fixed`,
with exactly the original requests, counts and design. Compare identical fits,
existing probabilities and expanded hypothesis families separately; BH values
may change because previously withheld tests enter the family. QL is unchanged.
No cohort or gene is selected for the follow-up based on its significance.

Both original and repaired call inventories remain descriptive observations
within this one selected study. The repair is not a calibrated-FDR or production
promotion claim.
