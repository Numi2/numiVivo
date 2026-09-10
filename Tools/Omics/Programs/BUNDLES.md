# Standalone native H5AD program bundles

`singlecell-h5ad-programs` scores supplied expression programs without adding
per-cell arrays to the pseudobulk JSON report. It uses the same native resolution
and arithmetic as the resident and pseudobulk routes, with flat resident
cell-by-program arrays and separately fingerprinted binary output. There is no
dense cells-by-genes allocation.

```sh
numivivo singlecell-h5ad-programs original.h5ad \
  --plan program-plan.json --output new-program-bundle
numivivo singlecell-h5ad-programs-verify new-program-bundle
```

The plan contains `schemaVersion: 1`, the usual H5AD `mapping`, the existing
`programs` definitions/budgets, `normalizationTarget` (default 10000), and
`matchFeatureNames` (default false). No cell selection or PCA fitting is implied.
All source counts contribute to normalization, including genes outside programs.
The default remains exact feature-ID matching. Explicit `matchFeatureNames: true`
uses exact original feature names; a requested name matching multiple source
features is rejected. Unrequested duplicate names are retained, all original
feature IDs stay in metadata, and every resolved source index is recorded.
There is no alias, case, ortholog or outcome-based conversion. The plan binds the
chosen mode and definitions retain their original supplied fingerprint.

## Artifact contract

| File | Meaning |
| --- | --- |
| `original.h5ad` | Independently owned source snapshot; APFS clone where supported |
| `plan.json` | Exact source mapping, program definitions/budgets, matching mode and normalization |
| `metadata.json` | Original ordered sample, cell and feature identities |
| `model.json` | Method, original resolved feature indices and weights, missing members, coverage, update count and source passes |
| `scores.bin` | Every cell/program coordinate, including zeros; `u32 row, u32 program, f64 score`, little endian |
| `detected-members.bin` | Every cell/program coordinate; `u32 row, u32 program, u64 detected count`, little endian |
| `total-counts.bin` | Every cell at column zero; `u32 row, u32 zero, u64 total`, little endian |
| `receipt.json` | Exact source, plan, metadata, model, three array and implementation fingerprints |

All binary records are 16 bytes and appear in complete row-major order.
**A total of zero means the score is unavailable.** Its binary score is a finite
+0 placeholder and detection count is zero. Consumers must apply this total-count
mask; a nonempty cell with zero score remains measured zero. Original JSON
routes continue to use null for an empty cell.

Publication performs two canonical sparse source scans: cell QC, then program
accumulation. It uses the shared 1 MiB buffered record writer and atomically
renames a private staging directory only after success. Metadata, cell QC and
flat cell-by-program arrays remain resident; this is not a wholly constant-memory
algorithm. Existing explicit score/update and source bounds still apply. It does
not raise the 512 MiB pseudobulk-report limit.

Verification checks the canonical receipt/plan and every recorded payload hash,
then independently snapshots the source and repeats the native computation.
A result with a newly forged payload hash must still reproduce the original
computation. The bundle does not infer labels or calibrate biological activity.

## Qualification

The owner tests cover flat arithmetic, signed/zero/missing binary rows, supplied
coverage and budgets, explicit normalization, name ambiguity and retained IDs.
The [CLI checker](check_bundle.py) generates independent CSR/CSC/dense fixtures,
checks every output coordinate/value against NumPy, compares the existing JSON
route exactly, checks repeat/replay and rejects corrupted and rehashed results,
ambiguous requested names, insufficient budgets, unknown scope and invalid
normalization/coverage. Fixture checks establish implementation behavior only.

The final release passes **12 tests in two suites** and **25 CLI checks**, with
13 expected rejections. The [qualification archive](evidence/2026-09-10-bundle/manifest.json)
retains all 234 stored/decoded members: fixtures, commands, outputs, build/test
logs, final modified source bytes, complete production hashes and the frozen
full-cohort protocol. The initial coordinator admission rejection and interrupted
collection attempt are retained. The latter copied no authoritative research
result; collection now explicitly excludes live private staging directories.

```sh
python Tools/Omics/Benchmarks/HIRISA/verify_archive.py \
  Tools/Omics/Programs/evidence/2026-09-10-bundle
```

The complete HIRISA protocol uses all 1,612,594 cells and both unchanged supplied
IFN definitions, with original Ensembl IDs and explicit exact-name matching.
Its [complete native publication, replay and comparison](../Benchmarks/HIRISA/NATIVE_PROGRAM_RESULTS.md)
against the previously frozen independent RNA reference now pass with exact
scores, detections, totals and metadata. This scale evidence is separate from
the smaller fixtures. The known integration-preservation failures remain
unchanged: these new bundles score original measured RNA, not corrected PCs.
