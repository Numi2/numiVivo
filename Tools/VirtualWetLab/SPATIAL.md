# Spatial tissue adapter v1

The same `WetLabExperimentAdapter` lifecycle now admits two assay families:
`cell-response` (RNA, unchanged verifier) and `spatial-tissue` (native spatial
transport). Each provides a catalog, prediction, reveal, summary and exact replay.
The registry rejects unknown families and record versions. Molecular perturbation
is a future adapter; RNA, tracer concentration and phenotype are distinct endpoints.

## What runs

NumiLab's `numi-matter-wet-lab` builds an authored tetrahedral scaffold using the
existing Matter material compiler and Metal FEM runtime. It applies a bounded
initial stretch, runs accepted mechanics steps and exports the accepted node
positions, source mesh digest and owner fingerprints. Rejected steps publish no
geometry. This is a one-way **frozen geometry snapshot**: mechanical time ends
before transport time begins. The example silicone constitutive law is a synthetic
scaffold material, not calibrated tissue.

NumiVivo's `numivivo-spatial` constructs extracellular volumes from those accepted
tetrahedra, subtracts explicit cell volumes, and places cell compartments at their
host tetrahedron centres. Shared faces carry `D × area / centroid-distance`
clearance; spherical cell surface area times declared permeability supplies
reversible passive membrane exchange. The existing physiological partition Metal
RK2 owner evolves concentrations with adaptive step reduction and conservation
certification. There is no Python prediction solver. This is a conservative graph
approximation; nonorthogonal continuum-mesh accuracy has not been qualified.

A plan fixes both initial fields, selected extracellular pulse regions, pulse
times, concentration increments and sample times. The matched control receives
no pulses. At each pulse boundary the native state is read, the declared amount
is added and a fresh native segment starts; the experiment retains global time
and the injected-amount ledger. Samples at pulse times are **after** the pulse.
External boundaries are closed. No advection, reaction, population growth, active
membrane biology, cell mechanics or feedback to Matter is implemented here.

The bundled example has 12 tetrahedra, 12 explicit cells, two pulses (0 and 2 s),
and samples at 0, 1, 2, 4 and 8 s. It is a numerical fixture, not a biological
specimen. Cell identity, position, host region and volume are explicit; a cell is
currently a well-mixed compartment, not a resolved cell surface.

## Build and register

Use absolute paths and new output directories. From the NumiLab checkout:

```sh
cmake -S matter -B /absolute/path/matter-build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build /absolute/path/matter-build --target numi-matter-wet-lab -j3
```

From NumiVivo:

```sh
python3 Tools/VirtualWetLab/build_spatial.py --output /absolute/path/vivo-build
```

The bounded build copies the exact partition owner Swift and Metal sources and
writes their SHA256 receipt. The full package also declares `numivivo-spatial`;
the acceptance used the bounded build, not a full MLX dependency build. Keep the
Swift resource bundle beside its executable and retain Matter's metallib at its
compiled path. Changing source, executable or declared shader fails exact replay.

```sh
python3 Tools/VirtualWetLab/prepare_spatial.py \
  --mesh Tools/VirtualWetLab/examples/spatial/mesh.json \
  --plan Tools/VirtualWetLab/examples/spatial/plan.json \
  --material /absolute/path/numi-lab/matter/materials/silicone.nmatter \
  --matter-binary /absolute/path/matter-build/numi-matter-wet-lab \
  --matter-shader /absolute/path/matter-build/shaders/NumiMatter.metallib \
  --spatial-binary /absolute/path/vivo-build/.build/release/numivivo-spatial \
  --partition-shader /absolute/path/vivo-build/.build/release/NumiVivoSpatial_NumiVivoShaders.bundle/NumiVivoPartitionKernels.metal \
  --output /absolute/path/spatial-assay
```

Optionally register `--observations observations.json`. A matching accepted
geometry must be generated first to bind the observation digest. The optional
`spatial_reference.py plan.json geometry.json observations.json` uses NumPy/SciPy
and an independent FP64 matrix exponential to generate **numerical-reference**
observations from the plan and geometry; it never reads native predictions.
Register these in a new assay directory, retaining the original registration.

```sh
printf '%s\n' '{"specimen":"two-block-scaffold"}' > selection.json
python3 Tools/VirtualWetLab/adapters.py predict /absolute/path/spatial-assay/assay.json \
  --workspace /absolute/path/experiments --selection selection.json
python3 Tools/VirtualWetLab/adapters.py reveal /absolute/path/experiments/RUN_ID
python3 Tools/VirtualWetLab/adapters.py verify /absolute/path/experiments/RUN_ID
```

NumiLab accepts repeated `--assay` flags, so this assay and Kang RNA share history
and the same interface. Spatial-only use does not require an RNA binary. Change
geometry, intervention fields or sampling by authoring a new registered plan;
the UI selects supported specimens/protocols instead of silently extrapolating.

## Record and observation contract

Registration precedes native execution. Each record copies mesh, material, plan,
accepted geometry and adapter sources. Matter must reconstruct the geometry before
NumiVivo predicts. The seal binds the registration, runtime identities, logs and
full two-arm prediction. Reveal never changes the prediction. Missing observations
remain unavailable; they are never replaced by the model's own output.

`wet-lab-spatial-observations/v1` requires specimen identity, geometry SHA256,
ordered compartment IDs, `units: mol/m3`, provenance, evidence class (`measured`
or `numerical-reference`), and both arms at every registered time. See the reference
writer for the complete JSON shape. Scoring retains RMSE and maximum absolute
error for every arm and time across **all** compartments. Wrong units, times,
identities, dimensions or nonfinite/negative fields fail. Even a measured-data
comparison does not automatically confer biological qualification.

Exact replay reconstructs both native owners and, after reveal, reconstructs the
comparison from bound observations. Hashes are integrity receipts, not signatures
or access-controlled escrow. Retain the complete record and exact native runtimes;
the UI's JSON download is only a summary. No random seed is needed for this
registered deterministic model; repeated execution is replay, not a biological
replicate. Each geometry/transport/cell/phenotype transition has its own evidence
status; this version does not represent a validated compound-to-tissue chain.

## Numerical acceptance

```sh
python -m unittest discover -s Tools/VirtualWetLab -p 'test_*.py' -v
NUMIVIVO_HDF5_LIBRARY=/absolute/path/libhdf5.dylib \
python Tools/VirtualWetLab/check_spatial.py \
  --assay /absolute/path/spatial-with-reference/assay.json \
  --output /absolute/path/new-acceptance \
  --rna-record /absolute/path/existing-rna-record --rna-binary /absolute/path/numivivo
```

[October 3 receipt](evidence/2026-10-03-spatial/checks.json) retains all nine cases:
baseline, half time step, disabled diffusion, disabled membrane exchange, inverted
geometry rejection, exact native replay, delayed reference reveal, tampering and
RNA compatibility. Maximum baseline concentration error was 3.416e-6 mol/m3;
half-step error was 7.954e-7. Baseline relative amount error was 1.181e-7; all
admitted variants were below 2e-5. These qualify this bounded graph calculation,
not continuum convergence, biological prediction or general mechanical coupling.

Development failures are retained: a Swift scoped-build entrypoint naming error,
a pre-existing partition pipeline lookup using a namespace despite explicit Metal
host names (fixed in its owner), and a test-journal exclusive-write error. The
journal was repaired and the suite rerun at unchanged model settings/thresholds.
