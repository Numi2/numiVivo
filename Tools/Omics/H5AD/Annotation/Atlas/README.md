# Full real HIRISA annotation

Native annotation now handles the retained 6.117 GB HIRISA atlas containing
1,612,594 observations. A source-bound plan adds `obs/numivivo_source_row`, with
each value equal to its original zero-based row position. This is ordering
provenance, not a cell-type, phenotype or prediction label.

The native run on an Apple M4 (10 CPU cores, 24 GiB RAM, macOS 26.6) completed in
6.26 seconds with maximum resident memory 189,677,568 bytes (180.9 MiB). This is
one measured local run, not a speedup comparison or cross-platform performance
qualification. Output logical size is 6,142,113,173 bytes versus 6,117,413,997
source bytes: growth of 24,699,176 bytes. The files have separate inodes.

## What changed

- Annotation admits 16 GiB inputs, 32 GiB outputs and up to two million rows or
  features. Existing total edit/payload limits remain in force; this does not
  admit an unbounded matrix or arbitrary numbers of new columns.
- Descriptor-rooted publication prefers APFS cloning, retaining separate content
  ownership and the existing atomic immutable/mutable publication rules. Other
  filesystems use bounded streaming. Forced streaming is independently tested.
- Detaching an edited HDF5 group copies its attributes and links unchanged
  children into the new group. Replacements change only the new group's links;
  aliases to original groups/datasets and old journal history remain intact.
- The shared raw-attribute storage allowance is set to annotation's output bound.
  The first atlas attempt stopped at its old 1 GiB default and published nothing;
  that failure remains in the evidence.

## Complete verification

A chunked independent HDF5 checker compared all 46 original stored datasets,
7,735,703,364 elements in total, including matrix data and structural indices.
Numeric bytes, string values, datatypes, shapes and original attributes agree.
Every added row value is correct. The complete original SHA256 remains unchanged,
and the output SHA256 agrees with the native receipt. The journal's canonical
plan hash and parsed plan agree with the submitted source-bound plan; the input
JSON file hash and native canonical-plan hash are explicitly distinguished.

The AnnData annotation suite passes on the final binary. Dedicated checks prove
hard-link group/dataset alias preservation, old journal retention and sharing of
unchanged child storage. The two-million-row boundary accepts and the next row
rejects; a sparse input just above 16 GiB rejects without publication. Two native
Swift tests and a standalone native check cover both preferred cloning and
forced streaming, byte preservation, independent inodes, mutation isolation,
permissions, immutable and mutable publication, symlinks, traversal, cancellation
and temporary cleanup.

The current scoped native build and all 106 recorded source inputs are retained.
This is not a full-product build, biological validation, new DE result, calibrated
annotation, Metal acceleration or evidence of general biological prediction.

## Reproduction and retained data

The [evidence manifest](evidence/2026-09-12/manifest.json) binds the source snapshot,
executables, plan, scripts, receipts, logs, failures and verification reports.
The original and annotated multi-gigabyte H5AD files remain locally retained at
the paths in `external-artifacts.json`, with sizes and SHA256 bindings. They are
not duplicated in the evidence archive. Their availability is required for full
reproduction; the archive alone is not a copy of the atlas.

With the archived scoped executable and plan, run:

```sh
h5ad-check annotate ORIGINAL_HIRISA.h5ad hirisa-plan.json NEW_OUTPUT.h5ad
python verify_atlas.py --source ORIGINAL_HIRISA.h5ad --output NEW_OUTPUT.h5ad --before hirisa-before.json --report verification.json
```

Use the recorded AnnData/h5py environment for the independent format checks.
The native annotation command requires no Python or scverse runtime.


## Full backed AnnData reopening

AnnData 0.13.3.post0 reopens the complete annotated file with `backed="r"` as
1,612,594 cells × 18,082 features, retaining a `_CSRDataset` backed matrix.
Both original indices and all 29 original observation/feature columns match the
source through the actual AnnData reader. Every new row value is correct.
Six deterministic count-row windows (14,931 stored entries) match exactly.
Count-window checking is sampled here; the earlier independent HDF5 check covered
all original stored elements.

The complete reference check took 26.85 seconds and peaked at 2,998,566,912 bytes
RSS. That includes metadata and reference-column comparisons; metadata remain
resident. It is not a comparable workload to native annotation and establishes
no native/scverse speedup. Both large file hashes were rechecked after the run.

[Backed-open evidence](evidence/2026-09-12-backed/manifest.json) retains the checker,
versions, full report and logs. Run `check_backed.py --source ORIGINAL_HIRISA.h5ad
--annotated ANNOTATED_HIRISA.h5ad --report backed-verification.json` in the recorded
AnnData environment. This closes the full-file reader check, not downstream
biological or prediction validation.
