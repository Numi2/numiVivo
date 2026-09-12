# Exposure and study identifiability

The executed design audit uses all 75 retained training donor identities and the
provisional study-level protocol durations: Kang 6 h, GSE181897 9 h, HIRISA 21 h.
It does not resolve the remaining donor/library exposure-admission gaps.

Three study indicator columns have rank 3. Adding duration produces four columns
but still rank 3. Duration is exactly reproduced as 6 times the Kang indicator plus
9 times the GSE181897 indicator plus 21 times the HIRISA indicator. The script checks
this equality exactly as well as computing numerical rank; the result does not
hinge on a floating-point rank threshold.

Thus an independent duration coefficient is not identifiable alongside unrestricted
study effects in this design. Removing study effects or imposing a mechanistic
response curve would be an explicit structural assumption, not new information
provided by the data. More donors at the same study-specific exposure do not add
within-study exposure variation. Dose/activity conversions have not been attempted.

`audit.py` accepts the retained model metadata and a new output path; `result.json`
records all donor rows, the null direction, ranks and source hash. No treated values
were read, model fit or predictions changed. This is an audit of the available
design, not experimental evidence for any causal explanation of failed transfer.
The next exposure-model data expansion should prioritize within-study exposure
variation and then independent evaluation rather than simply adding a study-constant
duration column to the failed model.
