# Held-out-animal count prediction with NB shrinkage

The frozen Crowell experiment finds a descriptive improvement in held-out
count prediction: the empirical-prior model gains **0.0426035 nats per gene**
over unshrunk MLE on the predeclared primary mean log score. Fixed-prior
shrinkage gains 0.0281420; empirical gains 0.0144615 over fixed. These are
plug-in NB scores conditional on observed sequencing depth, not calibrated
posterior predictive probabilities or proof of general effect-estimation benefit.

## Frozen design and source

[PROTOCOL.md](PROTOCOL.md) was frozen before fitting/scoring, SHA256
`07142ad0ead162b1bda03d78a847cca5b03a147084f430e6c15246eb4031ac87`.
The source is the complete deposited post-QC count scope already verified in
the [Crowell benchmark](../../Benchmarks/Crowell/README.md): 25,224 nuclei,
11,076 genes, 67,349,705 UMIs and eight independent mice (four per arm).
This retains the authors' upstream filtering and previously inspected cohort;
it is not a prospective untouched-study validation.

Every choice of one held-out Vehicle and one held-out LPS animal is evaluated
in every eligible population: **16 splits × 7 populations = 112 folds**.
Each training file contains six real pseudobulk rows and their original cell
membership identities, with three animals per arm. CPE cannot retain three
training animals per arm and remains unavailable. Metadata contain identities
and feature annotations; held-out counts are separate scoring inputs.

Training uses the existing native NB cohort, Gamma parametric dispersion trend,
count/presence/support gates and explicit library-size normalization. That
normalization choice was fixed before scoring and differs from the earlier
median-ratio comparison. The held-out size factor is its observed library total
divided by the geometric mean of training library totals. No held-out count
enters training normalization, dispersion fitting, gene filtering or prior
estimation. All three models share final training dispersions and refit nuisance
coefficients: MLE, a fixed normal contrast prior of one log2 unit, and the
empirical weighted-quantile prior.

Models are written before the independent native scoring command opens the
held-out count file. Scores compare the same training-admitted genes within a
fold. A failed MAP fit would make that fold's primary comparison unavailable;
the protocol does not remove inconvenient genes to obtain a passing comparison.
All 112 folds are available in this run.

## Animal-level results

Each actual animal's four held-out predictions are averaged first, followed by
equal weighting across eight animals and seven populations. The same animals
appear across populations and the training folds overlap. Genes, folds and
populations are not independent biological replicates; no independence-based
p-value or confidence interval is reported.

| Population | Training-tested genes across folds | Empirical minus MLE, nats/gene | Animals improved | Empirical minus fixed |
| --- | ---: | ---: | ---: | ---: |
| Astrocytes | 10,283–10,605 | +0.0325312 | 8/8 | +0.0057785 |
| Endothelial | 9,061–10,035 | +0.0325286 | 8/8 | −0.0023766 |
| Excitatory neurons | 11,051–11,066 | +0.0494019 | 7/8 | +0.0358356 |
| Inhibitory neurons | 10,852–10,966 | +0.0517888 | 8/8 | +0.0349264 |
| Microglia | 7,740–9,115 | +0.0362114 | 8/8 | +0.0004238 |
| OPC | 9,190–9,922 | +0.0440330 | 8/8 | +0.0104673 |
| Oligodendrocytes | 9,879–10,284 | +0.0517296 | 8/8 | +0.0161753 |

Empirical shrinkage improves the primary score in 55 of 56 animal/population
summaries. The exception is LPS animal LC026 in excitatory neurons (−0.000834).
Fixed shrinkage improves it in all 56. Empirical is worse than fixed in seven
endothelial, two excitatory-neuron and two microglial animal summaries; every
individual result is retained. These are 56 measurements of eight animals.

The secondary squared log1p normalized-count error decreases by 0.0113926 for
empirical versus MLE and by 0.0098647 for fixed versus MLE. Empirical has worse
secondary error than MLE in six animal/population summaries despite the overall
decrease. This predictive-score result does not establish effect-truth recovery,
posterior interval coverage, power, FDR control or a general automatic default.

## Native and independent verification

The measurement harness compiles the unchanged production Omics owner files at
source commit `e5ec73f80c69114553fbe249ee400aeab6dec7db`. It calls the existing
cohort, contrast-MAP and NB log-mass functions; there is no second estimator.
Its full-cohort probe exactly matches **31,950 product coefficient fits**,
including all MLE, fixed and empirical coefficients, SDs and scaled scores, plus
the design, dispersion trend and learned prior.

On the physical M4 Pro CPU, all **112 fits and 112 scoring commands** complete.
The 1,127,035 tested gene/fold pairs yield **3,381,105 saved conditional fits**.
Independent NumPy checks validate every training/test count against the source
aggregate, donor separation, training depth factors, feature filters, learned
prior, conditional score equations and information-based SDs. SciPy computes
every held-out NB log mass independently; the resulting aggregate native scores
agree to 3.55e-13 nats/gene. The largest scaled score is below 1e-7; the largest
SD difference is 6.72e-14 and secondary-error difference is 4.67e-15. There are
no final numerical failures or unavailable MAP fits.

The harness SHA256 is
`c7c0f9f0090685bc59ebece032a660c957dd1fba482c6400a459159ae831dab0`.
Four disjoint CPU processes completed the native fit/score work in 221 seconds.
This timing starts from prepared pseudobulks, excludes H5AD import and reference
checking, and is not a product end-to-end speedup or Metal measurement.

Retained preparation issues include a default-field JSON comparison failure and
a disk-full model mirror transfer, which stopped the first full reference pass
after 71 folds. Explicit default expansion fixed the former; streaming verified
models from the Mac mini removed the need for local mirrors. The complete
112-fold independent pass then succeeded. Original failed logs and attempted
checker sources remain. Verified regenerable Python bytecode and redundant,
hash-matched model mirrors were removed; complete native models and source data
remain retained.

## Reproduction

```bash
set -o pipefail
python Tools/Omics/NegativeBinomial/HeldOutRisk/prepare.py --root RUN --source_report CROWELL_FIXED_REPORT_GZ --product_report CROWELL_EMPIRICAL_REPORT_GZ
bash Tools/Omics/NegativeBinomial/HeldOutRisk/build.sh RUN/build
RUN/build/nb-heldout fit RUN/metadata.json.gz RUN/probe.json.gz | gzip -n > RUN/probe-model.json.gz
python Tools/Omics/NegativeBinomial/HeldOutRisk/check_probe.py --root RUN
python Tools/Omics/NegativeBinomial/HeldOutRisk/run_native.py --root RUN --binary RUN/build/nb-heldout --jobs 4
python Tools/Omics/NegativeBinomial/HeldOutRisk/check.py --root RUN
python Tools/Omics/NegativeBinomial/HeldOutRisk/summarize.py --root RUN
python Tools/Omics/NegativeBinomial/HeldOutRisk/archive.py --root RUN --out NEW_EVIDENCE
python Tools/Omics/NegativeBinomial/HeldOutRisk/archive.py --out NEW_EVIDENCE --verify
```

Reference dependencies are NumPy, SciPy and pandas. `check.py` also accepts
`--model-root REMOTE_RUN --model-host macmini` to stream the exact remote model
files. It checks them against the native output hashes before using them.
The full-cohort probe uses `nb-heldout fit metadata.json.gz probe.json.gz` and
the two source reports bound by the protocol.

[Evidence](evidence/2026-09-10/manifest.json) includes all training and test
inputs, native score outputs, independent checks, per-animal tables, plans,
source hashes and attempt logs. Complete fitted models remain at the external
paths in the manifest. This compact measurement format is not a new production
fit/prediction artifact API. Broader held-out-study risk, robust DE calibration,
posterior uncertainty and the full single-cell objective remain open.
