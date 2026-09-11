# Parse IFN-beta: full-cohort count validation

This is an input qualification for an independent expression-prediction test.
It does not yet establish predictive accuracy. The complete literal IFN-beta/PBS
cohort contains **725,031 cells, 12 donors, 24 donor/condition groups, 40,352 RNA
features and 1,373,870,697 stored count records**. No selected cell is downsampled
or removed because of annotation/QC disagreement.

**All twelve native donor ingestions and independent cell/aggregate count checks
now pass.** The retained matrix contains 3,070,817,047 counts, 32,500 fewer than
the historical source `tscp_count` sum. Historical detected-feature totals exceed
the retained matrix by 31,902; each historical QC field disagrees for 30,634 cells.
All 3,456 source runs are accounted for. [All twelve native source replays, final
payload retention and offline restoration now pass](COUNT_RESULTS.md). The full
result preserves every discrepancy and the earlier failed attempt.

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

The subsequent [HGNC identity audit](FEATURE_IDENTITY.md) finds 271 unique
approved/previous-symbol correspondence candidates among those 409 absences.
The other 138 remain unresolved or ambiguous, and the source lacks stable gene
IDs. These are nomenclature candidates, not a changed model panel or proof of
sequence-equivalent measurements; the exact-name compatibility gate remains open.

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
bash Tools/Omics/H5AD/build.sh "$NUMIVIVO_PARSE_STUDY/build-pooled" --with-cli
bash Tools/Omics/CountStore/Stream/test.sh "$NUMIVIVO_PARSE_STUDY/build-pooled"
python "$R/prepare_donors.py"
python "$R/test_transport.py"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$R/run_donors.py" --phase ingest
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$R/run_donors.py" --phase verify
python "$R/check_donors.py"
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
hashed and checked. This archive remains **preparation evidence**. The later [complete count result](COUNT_RESULTS.md)
separately retains the full count scan, independent comparison and native replay;
those results do not retroactively change the preparation archive.

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

## Memory correction, source failure and resumable execution

The original monolithic attempt stopped on an upstream HTTP 500 after
1,166,913,217 records (18,670,611,472 stream bytes; 2,396.62 seconds). The native
owner rejected truncation and published no complete bundle. That failure is
retained separately from the successful preparation and synthetic controls.

[Allocator controls](../../CountStore/Stream/Memory/README.md) identified temporary
memory growth and verified the per-chunk pooling correction. The source adapter
now retries only transport errors and HTTP 500/502/503/504, at most six attempts
with bounded backoff and the same `If-Match` range. Access and identity failures
remain fatal; three transport regression tests check that distinction.

The replacement run partitions the unchanged complete cohort by the twelve
original donors. Each partition retains every selected source cell and feature,
its original-to-local row map, native QC and two donor/condition aggregates.
Every native donor report is checked against independent count sums before its
completion receipt is saved. Resume validates completed artifacts and retries
only an unfinished donor. Full-cohort success requires all twelve ingest and
native verification receipts, exact disjoint coverage and matching source ranges.
This changes execution order, not cohort membership, biological replicates or
analytical endpoints. Each donor stream has its own SHA-256; a monolithic stream
SHA-256 is not claimed for the partitioned run.

For the current retained preparation, use the corrected `build-pooled` executable:

```sh
python "$R/prepare_donors.py"
python "$R/test_transport.py"
python "$R/run_donors.py" --phase ingest
python "$R/run_donors.py" --phase verify
python "$R/check_donors.py"
```

The first command creates a new immutable partition plan. Either execution
command can resume already completed donors. Native artifacts retain the exact
binary identity; changing the binary requires a separate qualification. The
full-cohort replacement ingestion and all twelve native replays have passed.
Neither is a predictive success or a completed biological validation.

The first donor completed native ingestion for 110,923 cells. The independent
checker initially decompressed an NPZ array once per cell; caching its arrays
fixed that validation bottleneck. Its saved native bundle and every independent
cell/aggregate count were checked without re-downloading. This recovery retains
the interrupted producer attempt and marks its unavailable peak RSS, exit status
and producer digest explicitly. Its separate full native source replay has now
passed, while the original missing process measurements remain unavailable.

## Retain the completed donor result

For this preserved qualification study and its frozen `donor-recipe-v2`, after
both phases and the complete-cohort check pass, retain the result with:

```sh
python "$R/retain_donors.py" --study "$NUMIVIVO_PARSE_STUDY" \
  --out "$R/evidence/2026-09-11-counts"
```

The retainer first requires all twelve ingestion and replay receipts. It checks
every retained preparation input against the immutable archive, reconstructs all
twelve donor plans from the original source axes and source ranges, validates the
executed recipe and binary identity, and reruns the independent full-cohort
checker. It writes a content-addressed archive, verifies every object and source
file again, and only then publishes the result directory. The exact existing
preparation and memory-control archives are required dependencies, preserving
axes, chunk maps and qualified executables without duplicating them. The dated
retainer requires the exact qualified executable already in its dependency
archive; a different build requires its own executable/source evidence and
separate qualification. Failed
source attempts and the first-donor checker recovery remain in the result.
The earlier incomplete cohort was checked to reject before creating any archive
directory. Full-result packing, archive verification and separate offline
restoration now pass; see [results and restoration commands](COUNT_RESULTS.md).

### Immutable preparation binding

The [binding evidence](evidence/2026-09-11-preparation-binding/manifest.json)
checks all 7,437 retained preparation inputs and 37 donor-plan files, including
all 3,456 cell-axis files and 3,456 count-range maps. Every donor plan must equal
the declared projection of the original parent metadata, with exactly the
original samples, features, cells, row cardinalities and source declaration.
Every source run and global-to-local row map is reconstructed independently
from the immutable parent; all 725,031 selected rows are covered exactly once.
Rehashing a substituted donor plan cannot make it an original source projection.

The complete check passes and repeats exactly. Seven substitutions against the
real Donor9 partition are rejected: barcode, condition, source declaration,
per-cell cardinalities, source-entry range, a rehashed row-map permutation and
total record count. Active inputs and the running binary were never altered.
The initial checker incorrectly counted historical range-receipt files as
canonical axes/maps; that failed check is retained. The corrected check requires
the exact canonical path sets and still hashes the historical receipts.

This qualifies input provenance, not count/replay or prediction completion.
The retainer rechecks all bound inputs after packing and publishes only if they
are still unchanged. To run the binding check independently:

```sh
python "$R/verify_donor_preparation.py" --study "$NUMIVIVO_PARSE_STUDY" \
  --preparation "$R/evidence/2026-09-11-preparation/preparation.tar.gz" \
  --out /path/to/new-preparation-bindings.json
python "$R/check_preparation_mutations.py" --study "$NUMIVIVO_PARSE_STUDY" \
  --out /path/to/new-mutation-check.json
```
