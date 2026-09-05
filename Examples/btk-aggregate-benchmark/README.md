# Published BTK occupancy aggregate comparison

`observations.json` records factual numeric summaries from Tam et al., *Phase 1 study of the selective BTK inhibitor zanubrutinib in B-cell malignancies and safety and efficacy evaluation in CLL*, Blood (2019), DOI **10.1182/blood.2019001160**. Primary record: https://pmc.ncbi.nlm.nih.gov/articles/PMC6742923/ . The relevant locator is the Pharmacodynamics paragraph describing paired lymph-node biopsies at week 1, day 3, predose.

The once-daily 320 mg cohort had 12 evaluable patients, median nodal occupancy 94%, range 82.4-100%, and a reported 50% above 95%. The twice-daily 160 mg cohort had 18, median 100%, range 86.3-100%, and a reported 89% above 95%. Percentages in the JSON are fractions. The rounded 89% is retained as reported; no exact patient count is reverse-invented from it.

These are published aggregate facts, not individual observations or recommendations. The source was read through indexed primary-article text; direct PMC access presented a browser challenge. No figure digitization, full-text redistribution, patient reconstruction, synthetic assay SD, or claimed PDF snapshot hash is included. The factual transcription is versioned by repository content identity; it does not establish authenticity of every future input or permission to redistribute the article.

## Inspect and compare

```sh
swift run numivivo occupancy-benchmark-inspect \
  Examples/btk-aggregate-benchmark/observations.json

swift run numivivo occupancy-benchmark-compare comparison-request.json \
  --output comparison-result.json
```

The second input is a `VivoAggregateBenchmarkRequest`: `benchmark` contains the published manifest and `predictions` supplies one matching `VivoPredictedOccupancyCohort` for every reported cohort. Each predicted cohort includes compound, target, tissue, visit, regimen label, explicit evidence for population modelling and timing mapping, and `cohortDraws`. Each row is one complete virtual cohort with the observed number of individuals. Posterior particles describing uncertainty for one individual are not a virtual population.

The evaluator reports expected cohort median, the distribution of cohort medians across draws, the fraction strictly above the reported threshold, and descriptive differences from the published aggregates. Published ranges are retained as ranges, not converted to confidence intervals or SDs. Missing cohorts are rejected, not omitted from a favorable score. There is no aggregate-to-individual conversion or Gaussian fitting likelihood.

The portable test uses explicitly synthetic all-one virtual cohorts to verify arithmetic and schema behavior. It is not a clinical prediction, not an actual fit to these data, and not an inference that either regimen meets a safety or efficacy objective.

## What remains missing

Individual occupancy measurements with assay precision, patient-level unbound intracellular exposure histories, a verified mapping from the clinical visit to simulated exposure/dose timing, and cohort-appropriate kinetic/population parameters are not supplied. The current code cannot produce an experimentally calibrated BTK prediction from this summary alone. Those gaps are recorded in the manifest rather than filled with undocumented estimates.
