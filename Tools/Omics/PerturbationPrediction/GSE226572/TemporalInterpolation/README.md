# Native donor-and-time-held-out response interpolation

A new time-aware baseline passes the predefined donor-level development gate on
GSE226572. It uses the same complete QC-admitted PBMC endpoint as the original
failed fixed-response transfer test: 126,633 cells aggregated into population RNA
profiles, 12,993 fixed shared features and all 18 treated donor/time observations.
It is not a B-cell endpoint or a new untouched external validation cohort.

For each observation, training excludes the entire query donor **and** every
observation at its query time in the other donors. The native Swift implementation
averages donor responses at each remaining time, interpolates linearly in hours,
uses a zero-response anchor at zero hours, and holds the final training response
constant after the final knot. All three 36-hour cases therefore use the final
24-hour response; no extrapolated slope is invented. Query controls are allowed;
query treated vectors never enter the native input. Nonnegative treated-expression
clipping follows response addition.

The time-independent baseline uses exactly the same training observations, weighting
donors equally and their available times equally within donor. No-change is the
query control. The protocol was written before executing this candidate, with no
parameter grid or outcome-selected time window. The gate requires at least 5%
lower equal-within-donor-time mean RMSE against both baselines in every donor.

| Held-out donor | Native RMSE | Training mean RMSE | No-change RMSE | Gain vs mean / no change |
| --- | ---: | ---: | ---: | ---: |
| D34 | 0.196497 | 0.234559 | 0.361566 | 16.23% / 45.65% |
| D38 | 0.202643 | 0.254525 | 0.401386 | 20.38% / 49.51% |
| D39 | 0.190859 | 0.246816 | 0.374268 | 22.67% / 49.00% |

All three donor gates pass. At individual times, the candidate is worse than the
training mean in **2/18** cases and worse than no change in **0/18**. Every case,
time knot and error remains in `results.json`; eighteen observations are not
eighteen independent donors. Three donors provide limited generalization evidence.
The original cross-study predictor's failed results remain unchanged and are not
superseded by this within-study model trained on different information.

The standalone native executable compiled and ran all 18 folds. All 233,874 native
candidate values and the corresponding baseline values match independent NumPy
interpolation/aggregation to at most 1.78e-15. Scalar compensated RMSE calculations
also pass. One complete native output repeats byte-for-byte. These checks are not
a full-product build, general temporal uncertainty calibration or biological
validation beyond this reused population-level endpoint. Cell composition and
preparation effects remain part of the endpoint.

Source arrays were checked byte-for-byte against the previously committed compressed
GSE226572 evidence. Only observed profiles and controls are used as training values;
previous model predictions are not training targets. Their SHA256 is
`534a0f5bcd676af37cda6af9b7085da012dc68edc90ac99c3fdffc06b971e341`.
`run.py` creates donor/time-excluded inputs, executes the native binary, freezes
output hashes before scoring, and verifies and scores all outputs. `evidence.tar.gz`
contains source, executable, all outputs, logs, protocol and manifest. Large repeated
input JSONs and the original source NPZ are external hash-bound dependencies.
Archive members passed independent size/SHA256 checks before publication.

Reproduce in a fresh workspace: compile `Temporal.swift` with
`swiftc -O -parse-as-library Temporal.swift -o temporal`, then execute `run.py` with
NumPy and one BLAS thread, adapting only workspace paths. Runtime evidence remains
at `/Users/n/numivivo-duration-interpolation-20260912`. No product prediction default
is changed. Next qualification needs additional donors and a distinct evaluation
cohort before promoting a general temporal predictor.
