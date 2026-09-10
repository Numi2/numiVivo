# Complete native HIRISA expression-program scoring

**PASS for full-source native publication, exact replay and independent numerical
comparison.** All **1,612,594 cells, 18,082 genes and 131 libraries** are retained.
Both unchanged supplied IFN definitions produce **3,225,188 cell/program scores**.
Every score, detected-member count and library total matches the previously
frozen [independent RNA reference](INTEGRATION_PROGRAMS.md); maximum score error
is **0.0**. Original metadata is byte-exact, preserving
Ensembl IDs, cell identities and source annotations.

This closes the native program artifact scale gap identified by the earlier
diagnostic. It does not change the integration-preservation failures or establish
pathway calibration, authoritative cell labels or biological outcome prediction.

## Scope and provenance

The [standalone native bundle](../../Programs/BUNDLES.md) uses two canonical
sparse source scans and shared program arithmetic. Metadata is separate from
binary scores, detections and totals. There is no dense cells-by-genes allocation
and no increase in the pseudobulk JSON limit. Metadata, cell QC and flat
cell-by-program arrays remain resident.

The release executable SHA256 is
`70a5bdcb258f35f947811bb7b5f2ddc777dfef19ce26880a74de818c8586ac95`.
It was built from `f1c1bc6a` plus the exact recorded native program patch, published
in `37e8f6e9`; the runtime environment retains every source-file hash. No later
commit identity is substituted for the actual build. Twelve native tests and
25 CLI controls, including thirteen expected rejections, preceded this full run.

The frozen source H5AD SHA256 is
`0873e698ebf8770a54e6dba09724ffbeda5e1a67dbf24d4e223552d4aca7969c`.
Exact original feature-name matching was explicitly selected without modifying
feature IDs or supplied definitions. All 93/97 matched Hallmark members and 5/5
five-gene members resolve to the same original indices and effective weights as
the independent reference. Four missing Hallmark members remain missing; the
five-gene set is overlapping evidence, not independent replication. Ambiguous
requested names are rejected. Unrequested duplicate names remain in the source.

All **3,845,991,249 nonzero entries and 7,700,096,227 UMIs** contribute to cell QC
and normalization. The run performs **79,599,318 member updates**. No cell has
zero total counts in this cohort. The separate fixture qualification verifies
that zero-count cells remain unavailable through the explicit total-count mask.

## Execution and verification

| Stage | Seconds | Peak resident bytes | Result |
| --- | ---: | ---: | --- |
| Native publication | 227.619 | 1982382080 | Complete |
| Independent all-cell check | 4.931 | 135872512 | Exact scores, detections, totals and metadata |
| Native source reconstruction/replay | 228.397 | 1976680448 | Exact receipt and derived payloads |

These are observations on the shared physical Mac mini, not controlled method
speed comparisons. The independent target arrays were already computed from
all original raw counts with a separate SciPy calculation; the short comparison
time does not include that earlier preparation. Its source, arrays, definitions
and execution freeze remain in the earlier integration-programs archive.

The native independent checker validates every explicit binary row/program
coordinate. It checks all resolved definitions/indices/weights, metadata/source
fingerprints, scores, detection counts and exact UInt64 totals, including
before/after artifact identities. Native replay independently snapshots the
original H5AD and repeats the two scans and publication computation.

## Reproduce and retain

The [native result archive](evidence/2026-09-10-native-programs/manifest.json)
contains all seven derived bundle files, including every score and original
metadata, plus protocol, exact commands, execution/source freezes, independent
checks, runtime environment and publication/replay logs. All 28
stored and decoded members are verified on both hosts. The original H5AD and
frozen executable remain externally preserved with exact hashes; independent
RNA arrays restore from the prior integration-programs archive.

```sh
python Tools/Omics/Benchmarks/HIRISA/verify_archive.py \
  Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-10-native-programs
```

`check_native_programs.py` checks a restored native bundle against the frozen
program reference and original metadata. `archive_native_programs.py` refuses
to archive a full result without completed native replay and independent checks.
Further annotation, rare-cell integration preservation and prospective
cross-context prediction remain separate requirements.
