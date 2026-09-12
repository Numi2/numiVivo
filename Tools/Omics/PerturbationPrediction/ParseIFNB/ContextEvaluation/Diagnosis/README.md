# Retrospective response-error diagnosis

The frozen candidate's applied response norm exceeds the observed norm in all
twelve donors (ratios 1.007–3.433). Response cosine is positive in eleven and
ranges from -0.217 to 0.640. Eight donors have greater error than no change.
The training-mean baseline has norm ratios 0.829–2.881 and cosine -0.247–0.711;
it is worse than no change in five donors.

Applied response means predicted treated expression after native nonnegative
clipping, minus query control. Excess MSE over no change is predicted-response
power minus twice its inner product with observed response, per feature.
Every decomposition matches direct vector and independent scalar sums within
1e-12. All twelve donors and all 11,600 features remain included.

This is post-hoc explanation of the retained failure, not a new validation test.
No scale was optimized, no donor was excluded and no prediction was changed.
Excess response magnitude coexists with imperfect direction; rescaling has not
been shown to solve transfer. Any calibration candidate must be selected using
training-only study holdouts and evaluated on untouched data. Parse outcomes
are now inspected and cannot supply a fresh independent validation claim.

[diagnosis.json](diagnosis.json) records every donor/model and source hashes.
[diagnose.py](diagnose.py) records the executed calculation and original host
paths. The original evaluation archive remains unchanged. This diagnostic was
created afterward, separately from its pre-score protocol.
