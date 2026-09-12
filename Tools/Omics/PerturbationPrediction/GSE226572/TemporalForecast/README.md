# Earlier-duration-only native prediction

This follow-up restricts the prior temporal baseline to information from strictly
shorter exposures in other donors. Every query donor is excluded completely; native
input validation rejects any training duration greater than or equal to the query
duration. The model applies the equal-donor mean response at the latest available
shorter duration, holding it constant beyond that knot. This is a simple persistence
baseline across exposure duration, not a fitted kinetic law.

All 18 original observations, three donors and 12,993 features remain. Two one-hour
queries have no eligible treated observation and explicitly return no change for
both model and training-mean baseline. They remain in scoring. The comparison mean
weights eligible donors equally and available times equally within donor. Clipping
and the whole-QC-admitted-PBMC endpoint remain unchanged from the previous experiment.

| Held-out donor | Candidate RMSE | Training mean RMSE | No-change RMSE | Gain vs mean / no change |
| --- | ---: | ---: | ---: | ---: |
| D34 | 0.235367 | 0.264120 | 0.361566 | 10.89% / 34.90% |
| D38 | 0.244206 | 0.277568 | 0.401386 | 12.02% / 39.16% |
| D39 | 0.228098 | 0.262503 | 0.374268 | 13.11% / 39.05% |

All three pass the pre-run requirement of at least 5% lower equal-within-donor-time
mean RMSE against both baselines. None of the 18 individual cases is worse than
either baseline; ties, including both fallbacks, are retained. The 233,874 native
candidate predictions and corresponding baseline values match NumPy to at most
1.78e-15, with independent scalar RMSE checks. All generated native inputs were
also inspected for strict donor/time exclusion. The standalone native build and
all 18 executions pass; this is not a full product build or default promotion.

The source experiment staggered stimulation so samples had a common culture and
harvest schedule. Shorter exposure therefore does not mean an earlier wall-clock
measurement. This test establishes only a duration-ordered information restriction,
not real-time or clinical forecasting. These three donors and outcomes were already
inspected; this is development reuse, not fresh biological generalization. Cell
composition, viability and preparation remain included in the population endpoint.
The original failed cross-study experiment is unchanged.

`Temporal.swift`, `run.py`, pre-run protocol, prediction freeze, all scores and
`evidence.tar.gz` retain the native experiment. The archive includes the executable,
all outputs and logs; each member passed size and SHA256 verification. Repeated
input JSONs and the original source NPZ remain external hash-bound dependencies
at `/Users/n/numivivo-duration-forecast-20260912` and the original retained source
path in the driver. The source is identical to the preceding interpolation test
and the committed original outcome artifact. Reproduce in a fresh directory by
compiling `swiftc -O -parse-as-library Temporal.swift -o temporal`, then running
`run.py` with NumPy and one BLAS thread. Additional donors and independent cohorts
remain required before a general temporal prediction claim.

The [all-feature audit](FeatureAudit/README.md) finds lower MSE than training mean
for 9,202 genes, but higher MSE for 3,769. Donor-average success must not be read
as uniform gene-level prediction accuracy. Every feature remains included.
