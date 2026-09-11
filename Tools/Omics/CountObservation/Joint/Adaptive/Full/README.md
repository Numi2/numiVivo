# Full-gene joint count fitting

The adaptive paired-rate model now runs across every source gene, with unavailable
and nonconverged cases preserved. This extends numerical fitting beyond the
16-gene development panel. It does not supply a new held-out biological outcome,
calibrated prediction intervals, or a general phenotype predictor.

## Completed Kang execution

All **15,706 source genes** and all **2,651 admitted training cells** were
processed. The source contains 1,479,543 retained nonzero counts across eight
paired donors. Full RNA depths remain the exposure denominator. Of 8,414 genes
with admitted control and treated cell dispersions, **8,405** reach the continuous
mean likelihood-gap target of `1e-6`. Eight hit `oracleLeafLimit`; one hits
`finiteSupportFitNotConverged`. The other **7,292 genes** remain unavailable under
the original cell-dispersion rules. None of these cases is admitted as a converged
prediction model.

The independent per-cell reference verifies all 15,706 records, **20,279,579
numerical values** and **3,986,925 final partition leaves**. The maximum scaled
discrepancy is `3.67e-12`. These checks cover count reconstruction, NB2/Poisson
likelihoods, rate derivatives, finite-support weights and moments, complete
partition coverage and the continuous upper bound. They do not measure biological
prediction error. The exact original and repaired runs, including every failure,
are retained in the [evidence manifest](evidence/2026-09-11-kang/manifest.json)
and [machine-readable results](evidence/2026-09-11-kang/results.json).

## Repair exposed by the full gene set

The baseline converged for 8,379 genes, with eight machine-resolution limits and
27 leaf-budget limits. Its support search waited for an evaluated witness above
`1 + tolerance`. When the witness was just below that threshold but its explicit
FP64 upper-bound allowance crossed it, subdivision could continue to machine
resolution without requesting another support point.

The actual Kang DNAJC8 counts reproduce that case: baseline lower score
`1.000000999592271`, upper `1.0000010015922718`, target `1.000001`. The repaired
search requests a support point once the witness exceeds
`1 + 0.9 * tolerance`, reserving headroom for the numerical allowance. Its final
continuous acceptance target remains `1e-6`; the finite-support tolerance remains
`1e-8`. DNAJC8 then converges with upper score `1.0000009569422676` and 516 leaves,
versus 3,650 leaves at the original resolution limit.

All **44 Swift tests in six suites** pass on the physical M4 Pro, including the
permanent [original-count DNAJC8 regression](DNAJC8-regression.json). This repair
changes the support-search trajectory; remaining failures and every before/after
status transition are retained. In particular, `symbol|PANK2` previously met
the bound but now exhausts the inner finite-support iteration budget. The repair
resolves 27 former limits and introduces this one new limit; it does not assert
universal optimizer convergence.

## HIRISA scope and remaining work

Preparation also retained all **18,082 HIRISA source genes**, **119,513 previously
admitted training cells**, and **291,097,122 nonzero counts** across five paired
donors. The source H5AD has 1,612,594 cells; the count model uses the previously
admitted cohort, not every source cell. Its full-gene native fitting and independent
verification are **still running** in four disjoint shard streams. No completed
HIRISA full-gene result is claimed in this release.

Training-only donor-exclusion needs cell dispersions re-estimated without the
excluded donor. Full-cohort dispersions cannot qualify a donor-held-out fit.
Mixing-distribution and dispersion estimation uncertainty, query-moment stability
beyond the original panel, and new independent treated outcomes remain open.
Historical HIRISA nominal 95% treated coverage remains **35.54%**. These training
fits did not score or use query-treated outcomes to change that result.

## Sparse execution and reproduction

`prepare.py` processes one donor/condition CSR group at a time, verifies its full
counts and depths against the published calibration, roundtrips CSR to CSC and
back exactly, and writes a compressed lossless CSC cache. It never builds a dense
cells-by-genes matrix. HIRISA preparation nevertheless peaks at **3,812,655,104
bytes**: this is a group-bounded transpose, not constant-memory ingestion.

`run.py` emits 64-gene NDJSON shards. A header supplies shared depth bins; each
record supplies only one gene's sparse counts. The native harness releases gene
intermediates after each result. Kang peaks at **19,906,560 bytes native RSS**
across its shards; this excludes Python preparation, reference verification and
the source cache. Per-shard receipts pin source/cache/manifest and
executable hashes. Every gene index is preserved, and interrupted partial files
require inspection before recovery. The immutable execution plan fixes shard size,
worker count and binary; worker locks prevent duplicate ownership.

The archive writer streams SHA-256-addressed files directly into parts of at most
48 MiB, avoiding a second whole-archive temporary copy. Verify or restore with:

```sh
python3 Tools/Omics/CountObservation/Joint/Adaptive/Full/retain.py verify \
  Tools/Omics/CountObservation/Joint/Adaptive/Full/evidence/2026-09-11-kang
python3 Tools/Omics/CountObservation/Joint/Adaptive/Full/retain.py restore \
  Tools/Omics/CountObservation/Joint/Adaptive/Full/evidence/2026-09-11-kang \
  /tmp/full-kang-evidence
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
NUMIVIVO_FULL_OUTPUT=/tmp/full-kang-evidence/study/Kang-gap-repaired \
python3 /tmp/full-kang-evidence/recipes/verify.py \
  /tmp/full-kang-evidence/study Kang
```

NumPy, SciPy and h5py are required. The retained cache is sufficient for the
independent reference and native replay. Rebuilding that cache from the original
H5AD uses the original calibration input metadata and native report; their paths
are recorded by `prepare.py`. Original source hashes and selected row identities
are retained in the manifest. The original H5AD is not duplicated in this archive.

Build into a fresh directory using `Tools/Omics/H5AD/build.sh`, then run
`Tools/Omics/CountObservation/Joint/Adaptive/test.sh` and `Full/build.sh` with that
directory. A decompressed `*-input.jsonl.gz` shard can be piped directly into the
resulting `full-counts` executable. For a new checkpointed remote run, create a
runtime JSON containing `fullBinarySHA256` and `nativePath`; set
`NUMIVIVO_FULL_RUNTIME` to that JSON and `NUMIVIVO_FULL_OUTPUT` to a new directory,
then invoke `run.py STUDY Kang`. Its SSH host is `macmini`. For four workers set
`NUMIVIVO_FULL_SHARD_COUNT=4` and distinct `NUMIVIVO_FULL_SHARD_INDEX=0..3` before
starting any worker. Do not change those settings for an existing execution plan.
`watch_verify.py` independently checks newly completed shards; `summarize.py`
requires complete source-gene coverage and a completed verification watcher.
