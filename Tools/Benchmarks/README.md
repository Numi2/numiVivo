# Independent molecular benchmark panel

The native `md-benchmark` command evaluates the existing Metal force kernels
against supplied observations. Python is used only to generate independent
OpenMM reference calculations, launch the CLI, and audit results.

```sh
python3 -m venv /absolute/reference-env
/absolute/reference-env/bin/python -m pip install -r Tools/Benchmarks/requirements.txt
/absolute/reference-env/bin/python Tools/Benchmarks/prepare_references.py --out /absolute/new-references
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
The vacuum protein–ligand complex retains its original vacuum model.

Limits and scientific boundaries are fixed in the
[campaign contract](../../Documentation/Design/FRONTIER_BENCHMARKS.md).
