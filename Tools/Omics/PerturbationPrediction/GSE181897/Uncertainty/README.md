# Why GSE181897 prediction intervals under-cover

**Observed cell sampling accounts for substantial variability that the original
training-donor intervals do not adapt to.** This diagnosis uses every admitted
donor and feature. It supports developing a count-based observation model with
explicit sampling information and zero-count uncertainty. It does not produce
a deployable interval or change the original external prediction result.

The [protocol](evidence/2026-09-11/protocol.json) was declared after all original
GSE181897 prediction scores were known and before this diagnostic. The work is
retrospective. In particular, the widened diagnostic below uses **observed
treated cells**, which are unavailable to the original prediction task.

## Actual cell resampling

All 4,938 admitted cells, 62 donor pairs and 11,800 frozen panel genes are included.
For each of the 124 donor/condition strata, 256 multinomial bootstrap replicates
resample whole cells with replacement at the original sample size of 4–122 cells.
Each replicate recomputes RNA normalization using all 20,303 source RNA features
before projecting to the panel. The seed and independent donor/condition streams
are fixed in the protocol. No cells, donors or genes are selected by the result.

The original H5AD aggregates match all 124 previously qualified count vectors.
For the first two replicates of every stratum, a separate repeated-row sum
matches weighted sparse multiplication exactly over every RNA feature: 248
independent count checks. A second complete run reproduces all moment arrays
and four result tables byte-for-byte. The two runs take 6.67 and 6.72 seconds;
these are diagnostic timings, not a controlled native performance benchmark.

## Coverage diagnosis

The original half-width is the frozen Student-t critical value times
`sqrt(trainingResponseVariance * (1 + 1 / trainingDonors))`. The diagnostic
adds observed control-cell and treated-cell bootstrap variances inside that
square root, preserving the mean, critical value, clipping and availability mask.
It is deliberately outcome-informed and has no calibration claim. Training
response variance already contains the measurement variation of its own libraries;
this addition is not a fitted decomposition of latent and measurement variance.

| Training origin | Original raw-response coverage | Widened raw-response coverage | Widened treated-expression coverage | Widened mean raw width |
| --- | ---: | ---: | ---: | ---: |
| Kang | 92.39% | 98.09% | 98.09% | 6.9664 |
| HIRISA | 31.01% | 89.72% | 94.25% | 3.9428 |

Widths use natural-log(1+CPM) units. The original HIRISA treated-expression
coverage is 35.54%; clipping accounts for the difference from raw-response
coverage. Its original mean raw width is 0.4282. Increasing width with access to
held-out outcomes is not a successful new prediction. Even this diagnostic
does not attain 95% raw-response coverage for HIRISA.

Re-estimating bootstrap variance separately from each half of the replicates
changes raw diagnostic coverage by at most 0.93 percentage points for a HIRISA
donor and 0.25 points for Kang. This checks finite simulation sensitivity; it
does not establish bootstrap convergence or biological calibration.

The minimum sampled cell count across a donor's control/treated pair correlates
with its raw-response RMSE at Spearman −0.945 for HIRISA and −0.951 for Kang.
These are descriptive associations, not causal estimates of sampling effects.

## Zero counts and residual structure

HIRISA raw-response coverage separates sharply by observed count support:

| Observed donor/gene stratum | Original coverage | Widened diagnostic coverage |
| --- | ---: | ---: |
| Zero in control and treated | 77.70% | 77.70% |
| Zero in control only | 0.87% | 64.67% |
| Zero in treated only | 0.32% | 78.43% |
| Nonzero in both | 22.09% | 99.47% |

All strata and both origins are retained in
[zero-strata.json](evidence/2026-09-11/results/zero-strata.json). Coverage uses the
original available-feature masks; these donor/gene pairs are not independent
biological replicates. Resampling observed cells assigns zero variance to a
gene absent from every sampled cell. It cannot represent unobserved transcripts
or cell states, so remaining undercoverage does not rule out sampling effects.

For the unclipped response residual `E = observed response - training mean`,
the exact squared-error decomposition is
`mean(E²) = mean_gene(mean_donor(E)²) + mean_gene(var_donor(E))`.
Here donor variance uses divisor 62, so the identity is exact.
The donor-average component is 4.74% of HIRISA's MSE and 15.95% of Kang's;
the remaining variation includes sampling, biological differences and other
effects. The mean observed-cell variance is 0.7663, about 58.62% and 51.72% of
the respective MSE magnitudes. These ratios are **not fractions of error causally
explained by sampling**. This decomposition precedes clipping and averages
squared errors, so it is distinct from the original equal-donor mean RMSE score.

## Next implementation requirement

Do not promote this outcome-informed widening or merely multiply the original
intervals until these known scores approach 95%. A usable extension needs:

1. A declared endpoint distinguishing latent population-average RNA from an
   observed finite-cell count profile, and explicit units for its uncertainty.
2. A source-bound observation model using training counts, available control
   cells and sampling information known at query time. Planned treated sampling
   depth or cell count must be explicit; donor aggregate rows must not be mistaken
   for individual source cells. Observed treated counts from this diagnostic must
   never become predictor inputs.
3. Count likelihood support for uncertain expression after zero observations,
   along with the distinction between latent donor variation and measurement
   variability. An arbitrary epsilon or observed-cell bootstrap alone does not
   resolve that support problem.
4. Validation of coverage, width, missingness and point estimates on new biological
   observations, with all unsupported contexts and earlier failures retained.
   GSE181897 is now development data for any such extension.

The [native count observation component](../../../CountObservation/README.md)
now implements the conditional NB2 likelihood and proper-prior posterior. It
preserves per-cell depths and distinguishes latent-rate uncertainty from future
sampling variance. Learned prior/dispersion, coupling to donor response and
independent calibration remain required before the original intervals are repaired.

## Evidence and reproduction

The [manifest](evidence/2026-09-11/manifest.json) binds all diagnostic arrays,
frozen model moments, original and repeated execution records, source identities
and environment versions. The actual recipe is [diagnose.py](diagnose.py).
It expects the existing qualified parent study at
`/Users/home/numivivo-gse181897-20260911`; rerunning the bootstrap requires that
source H5AD and its hash-bound count/role artifacts. The retained derived results
can be checked without those external source files:

```sh
python verify.py evidence/2026-09-11
```

NumPy, SciPy and AnnData/h5py versions are recorded in the evidence. For a fresh
full run, copy `diagnose.py` and the archived `protocol.json` into a new study
directory, preserve the parent source study, and run the script with the original
scientific Python environment. It refuses an existing `results` directory.
Native predictor code, existing fitted models and original predictions are unchanged.
