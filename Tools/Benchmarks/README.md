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
in checkpoints. Introduced in v5 and retained with the v6 RATTLE correction,
this opt-in mode currently supports fixed-cell classical
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

For smooth-Hamiltonian conservation tests, prepare a separately named reference
panel, then pass that new directory to the NVE runner:

```sh
/absolute/reference-env/bin/python Tools/Benchmarks/prepare_smoothed_references.py --references /absolute/new-references --out /absolute/smooth-references
python3 Tools/Benchmarks/run_nve_refinement.py --references /absolute/smooth-references --binary /absolute/frozen/numivivo --out /absolute/smooth-nve-refinement
```

Periodic LJ interactions explicitly
switch from 0.7 to 0.8 nm in both engines; independent observations are recalculated.
The vacuum case retains its original model. The original sharp-cutoff panel and
failed preparation remain available and must not be relabeled as passing this
different model. `nve_policy.json` and its limits remain unchanged.

The optional `plot_nve.py` helper uses `matplotlib==3.10.7` to plot the measured
energy deviations while showing missing cases explicitly. Install it in a
separate plotting environment, then run:

```sh
/absolute/plot-env/bin/python Tools/Benchmarks/plot_nve.py --qualification /absolute/smooth-nve-refinement/qualification.json --out /absolute/new-refinement.svg
```

Limits and scientific boundaries are fixed in the
[campaign contract](../../Documentation/Design/FRONTIER_BENCHMARKS.md).

For the longer conservation gate, pass `--policy Tools/Benchmarks/long_nve_policy.json`
to `run_nve_refinement.py`. It fixes 5 ps at each of three time steps and retains
the original conservation limits.

The separate interacting NVT campaign requires the Python environment above,
including the pinned NumPy/SciPy versions, and a frozen native executable with
its complete `NumiVivo_NumiVivoShaders.bundle`. Its JSON binary manifest contains
`buildSource`, `buildCommand`, and a `files` array with relative `path` and
`sha256` for the executable and every bundle file. No files may be omitted.

```sh
/absolute/reference-env/bin/python -m unittest discover -s Tools/Benchmarks -p 'test_ensemble_*.py' -v
/absolute/reference-env/bin/python Tools/Benchmarks/run_ensemble.py --references /absolute/smooth-references --binary /absolute/frozen/numivivo --binary-manifest /absolute/binary-manifest.json --out /absolute/new-nvt-campaign
/absolute/reference-env/bin/python Tools/Benchmarks/analyze_ensemble.py --campaign /absolute/new-nvt-campaign --out /absolute/new-nvt-qualification.json
```

The default policy fixes 12 rigid-water trajectories, independent seeds across
temperatures, two time steps, equilibration exclusion and block-bootstrap tests.
Analysis requires the recorded code and Python/NumPy/SciPy versions. It checks
the complete matrix, model identity, input/report hashes, clocks, independently
counted degrees of freedom and checkpoint constraints before using observations.
Missing data fail validation; insufficient statistical evidence stays inconclusive.
Parent preparation failures and excluded cases remain explicit. Therefore a
requested subset can pass while the full retained panel still returns nonzero.
The [ensemble design](../../Documentation/Design/MD_ENSEMBLE_QUALIFICATION.md)
defines the limited scientific claim.

An additional offline endpoint audit recomputes final potential energy with the
serialized OpenMM Reference system and start/end kinetic energy from saved masses
and exact velocity words. It retains the original per-particle potential limit;
the native v6 kinetic reduction is checked against its FP32 forward-roundoff bound.

```sh
/absolute/reference-env/bin/python Tools/Benchmarks/audit_endpoints.py --references /absolute/variant --campaign /absolute/native-campaign --out /absolute/new-endpoint-audit.json
```

Apply it separately to each NVE timestep or NVT matrix cell. This checks saved
endpoint accounting, not an independently integrated trajectory or experimental
agreement.

Native checkpoints wrap atoms separately. Before OpenMM evaluation, the auditor
reconstructs whole finite molecules using integer box translations and validates
every molecular cycle and exception image. Parameters, relative bonded geometry,
velocities and source checkpoints remain unchanged. Winding molecules, ambiguous
half-box bonds and unsupported cells are refused. This follows
[OpenMM's periodic-coordinate requirements](https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#periodic).

## Continuing an insufficiently sampled NVT campaign

The separate adaptive policy in `ensemble_continuation_policy.json` extends every
finer-step replica from its existing exact checkpoint to 210 ps, in 10 ps chunks.
It retains the first campaign's failures and inherits its scientific limits.
Run only after the original matrix and qualification have completed, using the
same executable/shader bundle and the same NumPy/SciPy/OpenMM environment:

```sh
python Tools/Benchmarks/continue_ensemble.py run --parent /path/to/nvt-campaign \
  --parent-qualification /path/to/nvt-qualification.json \
  --binary /path/to/frozen/numivivo --out /path/to/new-continuation
python Tools/Benchmarks/continue_ensemble.py analyze --campaign /path/to/new-continuation \
  --out /path/to/new-continuation-qualification.json
```

Output paths must be new. A failed chunk remains retained and ends that replica;
the other prescribed replicas are still attempted. Native execution and endpoint
agreement are distinct from the final statistical verdict. Source and policy
declaration alone do not establish a measured continuation pass.
