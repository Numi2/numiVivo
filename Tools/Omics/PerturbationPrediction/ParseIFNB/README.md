# Parse IFN-beta: full-cohort count validation

This is an input qualification for an independent expression-prediction test.
It does not yet establish predictive accuracy. The complete literal IFN-beta/PBS
cohort contains **725,031 cells, 12 donors, 24 donor/condition groups, 40,352 RNA
features and 1,373,870,697 stored count records**. No selected cell is downsampled
or removed because of annotation/QC disagreement.

The [protocol](PROTOCOL.md) fixes the source version, selection, checks and
scientific boundaries. The original H5AD is 227,497,986,816 bytes; an in-memory
metadata cache and bounded HTTP count ranges avoid retaining that full file.
Every response is version-pinned and hashed. The native [canonical stream owner](../../CountStore/Stream/README.md)
checks the count records and computes cell QC and biological-replicate aggregates.
This route retains cell metadata in memory; it is not fully out-of-core metadata
or a Metal performance claim.

The original `gene_count` annotation differs from retained X cardinality for
30,634 selected cells. Historical source QC is retained separately from matrix
QC; agreement is not forced by dropping cells or changing annotations. Missing
mitochondrial feature annotation leaves mitochondrial fractions unavailable.

A separate identifier check finds 12,584 of the duration model's fixed 12,993
response symbols after adding its documented `symbol|` namespace wrapper.
The other **409 symbols are absent by exact name**. No gene aliases are guessed,
missing features filled with zero, or panel silently changed. A future transfer
protocol must resolve correspondence or declare its supported panel before fit
and scoring; the current frozen model cannot simply be passed these axes.

The source methods identify 24-hour stimulation. The exact IFN-beta dose and
reagent remain unresolved in accessible primary evidence. No model has been
fitted or scored on this cohort. Count agreement would not repair the previous
external fixed-response failure, validate the duration model independently, or
establish predictions of tissue function or clinical outcomes.

## Reproduce

Use Python with NumPy and h5py, an Apple Swift toolchain, adequate RAM for the
metadata cache and cell axes, and internet access to the pinned source. Set an
absolute study directory. The cached index helper is specific to this exact
HDF5 version: it prefetches documented v1 index nodes, then uses HDF5's own
`chunk_iter` to enumerate and validate the entire chunk map. It never infers
count payload locations from file-size arithmetic.

```sh
export NUMIVIVO_PARSE_STUDY=/absolute/path/to/new-study
R=Tools/Omics/PerturbationPrediction/ParseIFNB
python "$R/bootstrap.py"
python "$R/read_roster.py"
python "$R/prepare_axes_cached.py"
python "$R/map_index_cached.py"
python "$R/assemble_plan.py"
bash Tools/Omics/H5AD/build.sh "$NUMIVIVO_PARSE_STUDY/build" --with-cli
bash Tools/Omics/CountStore/Stream/test.sh "$NUMIVIVO_PARSE_STUDY/build"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$R/run_counts.py"
python "$R/check_counts.py"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$R/run_counts.py" --verify
python "$R/check_counts.py"
```

The producer retains per-run range hashes and independent count totals, checks
its own completion and the native exit code, and bounds prefetch to sixteen
source runs. Replay checks every fetched range against the first execution and
requires native reconstruction. It does not retain the 21.98 GB canonical stream;
replaying requires the pinned source to remain available. An entire-file source
SHA-256 and validation of unselected counts are not claimed.

Original data: **Parse Biosciences, CC BY-NC 4.0**. Source extracts and derived
counts retain that attribution/license. Software is covered separately by the
repository license. See the [primary-source references and unresolved dose](PROTOCOL.md).

## Retained preparation evidence

The [preparation archive](evidence/2026-09-11-preparation/manifest.json) preserves
all selected axes, the complete HDF5-enumerated chunk maps, requested metadata
range identities, exact native executable, four passing native regression tests,
and failed/abandoned preparation attempts. All 103 scoped build inputs were
hashed and checked. The archive is **preparation evidence**: the full count scan,
independent aggregate comparison and native replay are still running/pending.
No complete count-validation result is claimed at this stage.

Verify or restore the archive with the existing exact-content archive reader:

```sh
python Tools/Omics/PerturbationPrediction/Duration/archive.py verify \
  --archive "$R/evidence/2026-09-11-preparation/preparation.tar.gz"
python Tools/Omics/PerturbationPrediction/Duration/archive.py restore \
  --archive "$R/evidence/2026-09-11-preparation/preparation.tar.gz" \
  --out /absolute/path/to/restored-preparation
```

The initial small-range metadata and repeated chunk-coordinate lookup attempts
were stopped because of request overhead. The completed route cached 1.60 GB
of tail metadata in RAM, then prefetched 14,724,400 bytes in 7,025 index nodes.
HDF5 enumerated all 393,213 data/index chunks; all 192 earlier per-run coordinate
maps matched. Index preparation took 160.79 seconds on this host/network; this
is not a matrix-ingestion or end-to-end performance result.
