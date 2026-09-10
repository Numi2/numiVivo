# Post-score dispersion-stage audit

Declared after the untreated-cell null benchmark at 6e6f8b2, before this audit's
stage-level reference fits. This is a diagnostic of that observed result, not
a new independent calibration experiment or a search for a passing filter.

Use all twenty original Kang/Hagai sham splits, their unchanged eligible genes,
paired designs and fixed native size factors. Refit pinned DESeq2 1.52.0 with
the same parametric/Wald/no-replacement settings as the frozen comparison.
Check every reproduced final dispersion, coefficient, standard error and
probability against the original reference table before interpreting any stages.
Record warnings, messages and failures; do not omit a split after failure.

Extract gene-wise and fitted dispersions, prior and residual log variances,
outlier admission and final dispersion for every eligible gene. Read all native
default reports through their existing compressed and logical SHA256 bindings.
Compare the original whole tested families and report the previously observed
twenty native Hagai calls separately. Keep native rank-withheld genes visible.

Use the stage comparison to distinguish gene-wise estimation, trend fitting,
prior variance and outlier admission. Do not change the native model, floor,
support policy, influence threshold or BH family in this audit. Agreement with
DESeq2 does not itself establish calibration, power or biological validity.
