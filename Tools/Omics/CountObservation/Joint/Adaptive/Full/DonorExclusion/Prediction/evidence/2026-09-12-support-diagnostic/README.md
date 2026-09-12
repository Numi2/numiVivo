# Frozen prediction support diagnostic

All 13 completed development folds were examined after scoring. The original scoring plan, predictions and primary results remain unchanged.

| Fold | Genes outside training rate range | Inside joint / baseline RMSE | Outside joint / baseline RMSE |
| --- | ---: | ---: | ---: |
| HIRISA-00 | 37.07% | 1.444 | 1.708 |
| HIRISA-01 | 30.38% | 1.770 | 1.704 |
| HIRISA-02 | 33.45% | 1.555 | 1.915 |
| HIRISA-03 | 45.71% | 1.332 | 2.224 |
| HIRISA-04 | 49.53% | 1.369 | 2.250 |

The joint model loses to the training-mean response baseline both inside and outside the range in all five HIRISA donors. Kang beats that baseline in both strata in all eight donors. The descriptive range mismatch alone therefore does not isolate the HIRISA failure. This is exploratory association, not proof of its cause or independent biological validation.

The range uses training donor negative-binomial MLE rates, whereas the query uses pseudobulk CPM; the estimators differ. No subgroup is promoted as a replacement benchmark. A revised model must be developed using training data only, frozen, and evaluated on fresh independent data.

The script checks the compressed score hashes, excludes the query donor from training IDs, and checks reconstructed control log1p CPM against frozen score records. JSON binds every inspected input shard by SHA256. Raw inputs and per-gene scores remain in the original study; this directory is a report and script snapshot, not a complete raw-data archive. To rerun, place the script in the prediction study root containing the 13 fold directories.
