# Fitting the count observation model from original training cells

`VivoCountObservationCalibration` learns gene-specific cell dispersion and a
proper Gamma prior for a new donor's latent RNA rate from source-bound training
cells. It connects those parameters to the existing
[count observation posterior](../README.md), replacing the fixed sensitivity
parameters used in that component's first numerical qualification.

Calibration keeps each donor and condition separate. It consumes sparse cells
through `VivoCellCountMomentAccumulator`, retaining arrays of sufficient moments
rather than a dense cells-by-genes matrix. Every cell includes its complete RNA
counts, so the native owner computes its full-library denominator before any
feature selection. This route supports training groups larger than the existing
512-observation donor-design GLM limit; it does not change that GLM's contract.

The method is a conditional moment estimator, **not** Cox-Reid likelihood fitting,
a calibrated Bayesian hierarchy or a repaired perturbation-response model.
Training parameter uncertainty is not integrated. A supplied condition identity
and distinct donor ID are required when applying the fitted model. Identifier
checks do not establish participant independence across studies.

## Within-donor cell dispersion

Let `e_i = fullRNA_i / 1e6`, `z_i = geneCount_i / e_i`, and `N` be the actual
number of observed cells in a donor-condition stratum. The model assumes
independent NB2 cells with common latent rate `r` and cell dispersion `phi`:
`Var(Y_i | r) = e_i*r + phi*(e_i*r)^2`.

For each gene and stratum, the native accumulator retains:

- Sample mean and variance of `z_i`, including all implicit sparse zeros.
- `P = mean(Y_i/e_i^2)`, an unbiased estimate of mean Poisson variance in CPM².
- `R2 = mean(z_i*z_j, i != j)`, an unbiased estimate of `r²` under independent cells.
- Exact counts, positive-cell count, full RNA total and observed cell count.

Thus `E[s² - P] = phi*r²`. Pool within-donor numerators and denominators with
weights `N-1`, then estimate `phi = max(0, sum((N-1)*(s²-P))/sum((N-1)*R2))`.
Donor mean shifts are not pooled into the cell dispersion. Negative raw estimates
are retained and labeled Poisson-boundary estimates. Zero denominators are
unavailable, including all-zero training genes and cases without sufficient
within-donor positive count pairs. Estimates outside the observation kernel's
supported dispersion range are unavailable rather than clipped.

The ratio estimator is not itself unbiased merely because its numerator and
denominator have unbiased expectations. Gene-specific estimates can be unstable
at low support; availability is an arithmetic admission result, not evidence of
accurate dispersion or calibrated intervals. No across-gene shrinkage is applied.

## New-donor rate prior

Each training donor contributes equally to the rate prior. The estimated
measurement variance of a donor's mean rate is `(P + phi*R2)/N`. Subtract the
mean of those variances from the observed between-donor sample variance, retaining
`tau² = max(0, observedVariance - meanMeasurementVariance)` and its boundary flag.
For `D` training donors, the new-donor rate variance is approximated by

`V = tau²*(1 + 1/D) + meanMeasurementVariance/D`.

This includes uncertainty in the estimated training mean. Even when `tau²=0`,
finite count information can provide positive variance; no arbitrary epsilon
is inserted. A Gamma prior with mean `m` and variance `V` has shape `m²/V` and
rate `m/V`. Nonpositive or unsupported parameter states remain unavailable with
a reason. The prior moment approximation conditions on estimated hyperparameters;
it is not an exact posterior predictive distribution for those hyperparameters.

The rate target is the average normalized rate across independent cells under
the fixed-depth NB2 model. It is not automatically identical to an RNA-weighted
pseudobulk CPM from heterogeneous cells. Rate/depth dependence, cell correlations,
mixtures and different assays can violate the model. The benchmark retains a
descriptive depth/rate association check without excluding genes by that result.
These distinctions must be resolved when coupling this layer to donor-response
prediction and defining the observed endpoint.

## Full original-data execution

The protocol admits every previously qualified training cell:

| Origin | Donors | Donor-condition groups | Original cells | Full RNA features |
| --- | ---: | ---: | ---: | ---: |
| Kang | 8 | 16 | 2,651 | 15,706 |
| HIRISA | 5 | 10 | 119,513 | 18,082 |

Both control and IFN-beta training conditions are fitted separately. No query
outcome is used. There is no gene filtering. Source feature correspondence,
cell barcodes, donor/condition membership and every complete RNA aggregate are
checked against the qualified parent cohort before accepting the stream.

The native owner runs on the physical Mac mini. The Python adapter reads original
H5AD CSR blocks, sends exact sparse cell records over SSH, and independently
accumulates reference moments. Cell records are streamed rather than retained as
a second multi-gigabyte matrix. The native reader rejects missing, duplicated or
unexpected source rows. Both sides hash every stream byte; the model fingerprint
binds the complete metadata bytes followed by the complete count stream.

HIRISA's `sampleID` is a preparation identifier, whereas the qualified cohort uses
GEO accessions. Admission explicitly verifies `geo_accession`, `geo_donor` and
`geo_treatment` for every selected HIRISA cell. An initial failed sample-column
assumption and the earlier nullable-string decoding failure are preserved.
The corrected adapter also supports AnnData nullable string indices without
accepting missing identities.

The reference reconstructs all moments and fitted parameters independently with
NumPy. Exact counts and support match exactly; numerical checks use
`abs(native-reference)/max(1,abs(reference))`, with limits 2e-8 for moments and
2e-7 for fitted parameters. All unavailable states and both boundary flags must
also match. Reference moment arithmetic differs from the native positive-cell
Welford accumulator and its direct distinct-cell product accumulation.

## Application to query controls

The fitted control models are applied to all 62 original GSE181897 control
donors and the 16 genes selected by identifier hash in the previous count
observation protocol: 992 requested gene/donor posteriors per training origin.
Unavailable calibration is preserved explicitly. Available cases use each
original control cell's RNA depth and counts. Planned sampling remains the
previously declared 20 same-condition cells with 1,000 RNA counts each.

This evaluates the complete native path from training cells, through fitted
observation parameters, to control-rate posteriors and future **same-condition**
count moments. It does not apply a treatment response or inspect treated query
counts. GSE181897 is known development data. Numerical agreement with an
independent posterior calculation does not establish improved biological coverage.

## Completed results

All **122,164 cells, 292,576,665 nonzero records and 553,598,348 RNA counts**
are processed. All 26 complete source aggregates match exactly. Independent
checks compare 1,728,464 moment values and all 67,576 gene/condition calibration
records. Maximum scaled errors are 8.65e-14 for moments and 2.64e-9 for fitted
parameters. All **29 regression tests in four suites** pass.

| Origin / condition | Available parameters | Insufficient count pairs | Dispersion outside kernel range |
| --- | ---: | ---: | ---: |
| Kang / control | 9,487 | 5,985 | 234 |
| Kang / IFNB | 9,262 | 6,230 | 214 |
| HIRISA / control | 14,709 | 994 | 2,379 |
| HIRISA / IFNB | 14,478 | 1,127 | 2,477 |

The range failures are computational admission limits, distinct from missing
count information. Raw estimates, including large dispersion and negative
Poisson-boundary values, are retained. They must be assessed rather than hidden
by clamping, default priors or reporting only successful genes.

For the 1,984 requested query cases, **1,178 posteriors are available and 806
remain unavailable**: 496 Kang cases lack within-donor count pairs; 310 HIRISA
cases have out-of-range dispersion. All available posteriors and their planned
count moments agree with independent adaptive quadrature; maximum relative
posterior error is 2.55e-12. All 371 available all-zero query cases retain positive
posterior variance. No treated outcome is used or predicted in this application.

An initial HIRISA run exposed excessive temporary-buffer retention in the
qualification reader: native peak RSS reached 4,701,716,480 bytes. An
`autoreleasepool` per cell limits those Foundation read-buffer lifetimes.
Complete replays preserve **every native report byte** and all independent
reference moments:

| Origin | Initial native peak RSS | Final native peak RSS | Initial streamed wall time | Final streamed wall time |
| --- | ---: | ---: | ---: | ---: |
| Kang | 266,665,984 B | 214,040,576 B | 15.27 s | 6.29 s |
| HIRISA | 4,701,716,480 B | 283,688,960 B | 686.04 s | 218.97 s |

HIRISA native peak RSS falls **93.97%**. The revised transport also enables SSH
compression. Wall times include transport and adapter work, and were not a
controlled speed benchmark; memory is the native process's measured RSS, not
whole-pipeline memory. The final HIRISA native report is 35,200,275 bytes; the
Kang report is 26,589,323 bytes. These remain complete, not reduced summaries.

Initial training executable hashes were not captured before a source-identical
runner relink. That provenance limitation is retained. The corrected training
binary is hashed **before and after both complete replays**, with full result
byte equality; the query binary is likewise recorded before application. The
final exact executables, source snapshots, both attempts, queries and reference
results are retained in the [evidence manifest](evidence/2026-09-11/manifest.json)
and [results](evidence/2026-09-11/results.json).

The prior remains specific to its source assay and population. The studies
normalize over 15,706, 18,082 and 20,303 measured RNA features respectively;
matching a gene symbol does not establish equivalent assay denominators or
valid cross-study prior transfer. Numerical qualification does not establish
biological coverage or repair the previous treatment-response intervals.

## Reproduction

Build the scoped runtime and numerical runners:

```sh
bash Tools/Omics/H5AD/build.sh /new/build --with-cli
bash Tools/Omics/CountObservation/Calibration/test.sh /new/build
```

The archived recipes preserve the original study paths. `prepare.py STUDY`
freezes cell identity and source bindings; `stream.py STUDY ORIGIN` reads the
original source counts and executes the native training tool on the Mac mini.
`verify.py STUDY ORIGIN` checks the complete output. `prepare_queries.py STUDY`
binds the fitted models to original control vectors; `Apply.swift` executes
those queries, and `prepare_reference.py STUDY` retains all unsupported cases
while preparing available posteriors for the independent reference in
`Tools/Omics/CountObservation/verify.py`.

The lossless archive stores each unique byte sequence once and maps all logical
paths to those objects, including identical complete replays. Restore and check
it in a new directory:

```sh
python Tools/Omics/CountObservation/Calibration/retain.py restore \
  Tools/Omics/CountObservation/Calibration/evidence/2026-09-11 /new/restored
python /new/restored/recipes/verify.py /new/restored/study Kang
python /new/restored/recipes/verify.py /new/restored/study HIRISA
```

All 85 logical files and 71 unique objects were restored and verified. Original
H5ADs and parent source-qualified cohorts remain hash-bound external inputs for
reconstructing the raw cell streams. Fitted models, every source moment, query
cell vectors, full results, exact final executables and both execution attempts
are retained in the archive. The initial source snapshot and final runner
recipes are separately retained so the memory repair is not lost.

The statistical APIs are part of NumiVivoKit. The framing and study adapters are
qualification tools, not new production CLI commands. Source files and decoded
metadata must remain immutable during the streamed execution; fingerprint and
membership checks are necessary for a source-bound result.

## Remaining prediction work

Integrate observation likelihoods with a latent donor-response model, preserving
the distinction between training measurement variance and donor variation.
Resolve the endpoint and depth assumptions before replacing the previous
pseudobulk intervals. Propagate or assess parameter estimation uncertainty,
measure behavior of unsupported genes and validate coverage and useful width
on new biological observations. The original GSE181897 treated-interval failure
and the rest of the twelve-part development goal remain open.

The subsequent [paired-donor diagnosis](../Paired/README.md) completes all-gene
and all-donor-omission checks. Independently corrected marginals do not guarantee
a valid joint covariance: 7,489 Kang and 6,373 HIRISA genes fail that condition.
The full endpoint comparison and covariance failures are retained before any
treatment-response connection is promoted.
