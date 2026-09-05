# Assay uncertainty, finite drug and measurement selection

All examples in this directory are synthetic mathematical inputs. They are not measured BTK kinetic parameters, prescribed clinical exposures, or treatment recommendations. The separate `Examples/btk-aggregate-benchmark` directory contains a small published-data summary with explicit missing-input boundaries.

## Portable checks

From the root of a complete checkout:

```sh
bash Tools/Posterior/run_portable_checks.sh /tmp/numivivo-research-checks
```

The retained Linux run passed 94 scoped checks. This compiles selected numerical code with a harness artifact adapter, not the complete package. See `Tools/Posterior/README.md` for exact compilation and verification scope.

## Application workflow on Apple

The main package must compile on the supported Apple toolchain. These complete-package command invocations were not run in the Linux environment; their underlying numerical operations were exercised by the portable harness.

```sh
swift run numivivo engagement-fit \
  Examples/target-engagement/synthetic-assay-inference.json \
  --checkpoint /tmp/assay-checkpoint.json --output /tmp/assay-posterior.json

swift run numivivo engagement-predict /tmp/assay-posterior.json \
  --output /tmp/assay-predictions.json

swift run numivivo engagement-design /tmp/assay-posterior.json \
  --candidates Examples/target-engagement/design-candidates.json \
  --output /tmp/measurement-design.json

swift run numivivo finite-drug-run Examples/target-engagement/finite-drug.json \
  --output /tmp/finite-drug-result.json

swift run numivivo occupancy-benchmark-inspect \
  Examples/btk-aggregate-benchmark/observations.json
```

Existing result files require explicit `--force` to replace. All operations can publish through the existing store with `--store <artifact-root>`. No CLI command silently switches to a GPU likelihood, supplies unreported assay uncertainty, or performs a laboratory experiment.

## Example meanings

`synthetic-assay-inference.json` uses two calibration conditions and one held-out condition. Generating association and dissociation rates are 100000 M^-1 s^-1 and 0.2 s^-1. Values come from the analytic reversible-binding solution plus fixed synthetic residuals, not from the NumiVivo forward routine under test. Reported measurement SD is deliberately absent. A third fitted parameter explicitly supplies residual SD with an assumed log-uniform prior from 0.002 to 0.04 fraction units. Absence of a reported SD is never relabelled as a measured precision.

Optional `assays` associate a case with a Gaussian noise model and a map of observation identifiers to exact/below/above/interval support. The fit fields `assayNoiseScale`, `assayNoiseFloor`, `assayBias`, `assayCorrelationFraction`, and `assayCorrelationTime` use the same explicit sharing and bounded-prior machinery as kinetic parameters. The noise model combines scaled reported variance and additional residual variance, plus an optional bias and exponential temporal correlation. Correlated censored observations are rejected rather than treated as independent. An interval describes one observation, not a cohort minimum/maximum.

`design-candidates.json` proposes a single measurement at each of four times under the explicitly selected same-context template. No outcome is supplied. Ranking estimates information about the whole joint fitted parameter vector, including assay nuisance parameters, with conditional Monte Carlo error and cost normalization. It does not claim a globally optimal experiment or design a multitime correlated panel.

`finite-drug.json` starts with 1e-6 M free drug and 2e-6 M free target and includes reversible/covalent binding, turnover, metabolism and clearance. Its evidence labels are explicitly assumed. The result is a local fixed-volume state transition with drug and target material ledgers. It is not a PBPK model, a dose schedule or a patient response. Transport and coupled commit remain the responsibility of the owning simulation.

The earlier mean-only `engagement-sensitivity` path intentionally rejects extended training assay models. A covariance/censoring-aware information calculation is needed before its rank could be used for these fits; the old J'J formula must not silently omit learned variance parameters.

Numerical uncertainty, assay uncertainty, population variation and missing mechanisms remain different quantities. This increment does not implement hierarchical donor inference, active-metabolite binding, a calibrated cellular endpoint, or GPU posterior evaluation.
