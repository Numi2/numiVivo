# Independent molecular benchmark panel

The native `md-benchmark` command evaluates the existing Metal force kernels
against supplied observations. Python is used only to generate independent
OpenMM reference calculations, launch the CLI, and audit results.

```sh
python3 -m venv /absolute/reference-env
/absolute/reference-env/bin/python -m pip install -r Tools/Benchmarks/requirements.txt
/absolute/reference-env/bin/python Tools/Benchmarks/prepare_references.py --out /absolute/new-references --project-initial-constraints --position-precision compensated
swift build -c release --build-tests --jobs 3 -Xswiftc -enable-testing
python3 Tools/Benchmarks/run_campaign.py --references /absolute/new-references --binary .build/release/numivivo --out /absolute/new-campaign
```

The directories must be new. Each reference directory contains original source
inputs, their provenance, serialized OpenMM System XML, and the complete native
request with three coordinate/force/energy observations. Each campaign case
retains its exact command, stdout, stderr, exit status and native report. The
scorecard independently recomputes comparison metrics and records failed,
unsupported and inconclusive evidence without converting them to successes.

The reference command returns nonzero if a source cannot be prepared; the
campaign still runs all successfully prepared cases and retains the missing
case. The original `water` AMBER restart has a skew cell incompatible with the
fixed 0.8 nm reference cutoff. `water-orthogonal` is a separately named OpenMM
package water fixture; it does not replace or repair the original source.

Use `--cases water-orthogonal alanine` for a small reference subset, or `--steps 0`
to request static comparisons only. The default 100 steps are a short execution
smoke check, not a thermostat/barostat equilibrium or performance qualification.
`--project-initial-constraints` explicitly prepares the imported geometry on
its declared constraint manifold before assigning velocities; the resulting
prepared checkpoint is retained. Static reference comparisons still use the
original supplied coordinates. Omitting the flag preserves the original state
and can correctly reject geometries that do not meet the execution constraints.
The vacuum protein–ligand complex retains its original vacuum model.

Compensated positions preserve GPU coordinate corrections and both exact words
in checkpoints. This opt-in v5 mode currently supports fixed-cell classical
dynamics with physical atoms. The original FP32 configuration remains available
for comparison; it fails the strict realistic-system constraint gate.

To reuse verified reference physics with a different native execution setting:

```sh
python3 Tools/Benchmarks/derive_campaign.py --references /absolute/new-references --out /absolute/variant --position-precision compensated
python3 Tools/Benchmarks/run_nve_refinement.py --references /absolute/new-references --binary /absolute/frozen/numivivo --out /absolute/nve-refinement
```

The derivation retains parent hashes, preparation failures, identical reference
geometries, models, observations and limits. The NVE runner copies and hashes
`nve_policy.json` before running three matched-duration timestep variants. It
checks identical prepared coordinates/velocities, energy conservation and
refinement using a bounded observation series. This is a short conservation
study, not equilibration or general scientific qualification. Timings with
observation collection include that work and are not throughput comparisons.

Limits and scientific boundaries are fixed in the
[campaign contract](../../Documentation/Design/FRONTIER_BENCHMARKS.md).
