# Target engagement: exposure, reversible binding and covalent conversion

This example connects an explicitly sourced kinetic model to unbound drug exposure and target occupancy over time. It includes a portable FP64 reference, lowering into NumiVivo's existing F1 ProgramPack compiler, an existing-runtime Metal cohort runner, a physiology adapter, and held-out observation evaluation.

**All supplied numbers and observations are synthetic mathematical fixtures. They are not measured BTK/zanubrutinib parameters, pharmacokinetic data, or treatment instructions.**

## Run the portable checks

From the repository root, with Swift 6 installed:

```sh
bash Tools/Kinetics/run_reference_checks.sh
```

This compiles five Foundation-only implementation files plus the check harness. It does not build the full package, link Apple frameworks, compile C++/Metal, or run a GPU. The recorded Linux run passed 44 assertions. See `Tools/Kinetics/portable-results.json` for exact source identities and the execution boundary.

Fixture regeneration is explicit and overwrites files in the selected directory:

```sh
bash Tools/Kinetics/run_reference_checks.sh --generate-fixtures /tmp/numivivo-kinetic-fixtures
```

## Run the integrated CLI on Apple silicon

The following commands require the main package to build on an Apple SDK. They are provided for Apple qualification; that package-wide build was not run during this development pass.

```sh
swift run numivivo engagement-validate Examples/target-engagement/synthetic-pulse.json

swift run numivivo engagement-run Examples/target-engagement/synthetic-pulse.json \
  --backend reference --output /tmp/engagement-reference.json

swift run numivivo engagement-source Examples/target-engagement/synthetic-pulse.json \
  --output /tmp/engagement-program.json

swift run numivivo engagement-compile Examples/target-engagement/synthetic-pulse.json \
  --output /tmp/engagement-program.nvpack

swift run numivivo engagement-run Examples/target-engagement/synthetic-pulse.json \
  --backend metal --output /tmp/engagement-metal.json

swift run numivivo engagement-study Examples/target-engagement/synthetic-study.json \
  --output /tmp/engagement-study.json
```

Existing outputs are rejected unless `--force` is supplied. `--store /path/to/artifacts` also publishes through the existing content-addressed artifact store. Backend selection is explicit: a Metal failure never silently becomes a CPU result. `engagement-batch` accepts a JSON array of experiments with shared kinetics, initial fractions, and sampling times but different exposure schedules.

The reference run record retains the experiment, numerical policy, result, and their fingerprints. The Metal batch result retains the actual program identity, numerical policy, device name, and experiments. Neither record is a clinical validation certificate.

## Conditional free-energy conversion

```sh
swift run numivivo engagement-rate Examples/target-engagement/synthetic-rate-request.json \
  --output /tmp/conditional-rate-evidence.json

swift run numivivo engagement-apply-rate Examples/target-engagement/synthetic-pulse.json \
  --rate Examples/target-engagement/synthetic-rate-request.json \
  --store /tmp/numivivo-artifacts --output /tmp/derived-rate-experiment.json
```

The conversion accepts an activation Gibbs free energy relative to a pre-reactive bound complex and an explicit classical transmission probability. It rejects an electronic-energy difference alone or an incompatible reference state. It does not calculate electron structure or sample a reaction pathway. Applying a derived rate changes only inactivation; association, dissociation and turnover remain unchanged. The derivation is stored before publishing the new experiment. Assumed inputs remain assumed, and total kinetic uncertainty remains unknown even when a conditional barrier uncertainty was supplied.

## Apple integration checks

```sh
swift run --package-path Examples/target-engagement target-engagement-checks
swift run --package-path Examples/target-engagement target-engagement-checks --metal
```

These check ProgramPack compilation, execution contracts, initial and recurring surrogate refresh, stale/failed receipts, actor single-flight behavior, missing Core ML uncertainty, bounded document I/O, and optionally Metal/reference occupancy agreement. They were syntax-parsed but not built or executed in the Linux development environment. A successful future run should retain its actual device/SDK/compiler information and output; the presence of this executable is not a recorded pass.

## Model boundary

Unbound ligand is an externally maintained concentration, not a finite pool depleted by target binding. Four normalized target states are supported: free, reversibly bound, covalently bound, and bound to one optional reversible competitor. All states share one turnover rate; target synthesis replenishes the baseline abundance. Exposure is right-continuous and piecewise constant, with no hidden interpolation or extrapolation.

The output is target occupancy, not cell survival, tissue response, toxicity, or an individualized treatment prediction. The study evaluator is not a parameter fitter or a posterior sampler. It rejects declared group leakage across partitions and retains failed cases instead of reporting favorable metrics over survivors.

For equations, source locations, compatibility changes, and integration details, see `Documentation/Design/TARGET_ENGAGEMENT.md` at the repository root.
