# Unseen-target co-response reference protocol

Frozen before scoring, 2026-09-09. Use the pinned full Norman reference from
`combinations.py`. For each of 105 targets, exclude every single and paired
condition involving it. Fit only control and remaining singles. Score the held
single; pairs and independent biological contexts remain future qualifications.

Descriptors are each target gene's log1p-CPM responses across available training
single perturbations. Use unique exact source-symbol matches only; preserve
unsupported aliases. The held gene's readout under other interventions is allowed;
its own intervention is excluded. Full genes remain in CPM denominators.

Fit multi-output ridge with intercept, alpha=1 fixed. Center descriptor coordinates
across training targets; divide by population standard deviation (constant: 1)
and sqrt(coordinate count). Predict all genes with dual ridge; compare with
scikit-learn SVD Ridge. A shuffled control rotates sorted training targets' output
responses by one. Self-intervention entries in training descriptors remain: their
absence for the held target is an acknowledged distribution mismatch.

Compare noChange, meanSingleResponse (all available singles), coResponseRidge,
and shuffledCoResponseRidge. Clip negative predicted log-expression to zero,
record clipping/implied CPM sums, and do not reclose. Persist frozen predictions
before separate scoring. All-genes and training/control mean-CPM top-1000 panels
(Ensembl-ID ties) are fixed without held outcomes. Report response RMSE/MAE/Pearson
and sign agreement; retain failures and matched-supported-target comparisons.

Check every fold's inputs after mutating all excluded conditions, independent
ridge agreement, and a complete repeated run. No outcome-based tuning or gene
activation filtering. This is pooled K562 CRISPRa transfer, not independent-donor,
Bayesian, mechanistic, or native-runtime qualification.
