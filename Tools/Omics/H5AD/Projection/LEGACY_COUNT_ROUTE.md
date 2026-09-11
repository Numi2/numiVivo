# Complete original Kang analytical input route

Native legacy projection, source-bound annotation, streaming count-store import
and pseudobulk aggregation now execute consecutively on the original Kang file.
All **24,673 cells × 15,706 genes** are retained. No Python-written replacement
H5AD or reduced cohort is used in this route.

The first attempt exposed a real handoff failure: the legacy projector preserved
`uns` without its dictionary encoding attributes, and native annotation rejected
that parent. Legacy projection v2 now adds the missing `dict/0.1.0` tag while
preserving its children. Existing modern input behavior is unchanged. The v1
source, binary, successful AnnData comparisons and subsequent failed annotation
attempt remain retained; they are not relabeled as v2 evidence.

## Inputs and interpretation

The original source SHA256 is
`e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830`.
Its float32 X values are all finite, nonnegative integers, ranging from 1 to
3,828. Recomputed total RNA counts and detected features match `nCount_RNA` and
`nFeature_RNA` for **every cell**. SCT totals match only 298 cells and SCT feature
counts match 13,403; those alternate assay annotations are preserved rather than
used to reinterpret X.

The explicit mapping retains the earlier source-qualified Kang UMI designation,
with `matrixPath=X`, eight original donor identities and the source `ctrl`/`stim`
conditions. The checker joins original `replicate` and `label` values with `__`
only to prepare a source-bound annotation plan. Native annotation adds that
`native_sample` column and a provenance journal; all original fields remain.
All sixteen sample definitions use the original donor as biological replicate.
Batch stays `unreported`. Author cell types remain analysis strata; this route
does not infer new labels or independent biological replicates.

## Verified results

| Check | Complete result |
| --- | --- |
| Native count records | All 14,184,532 row/feature/count triples agree with independently canonicalized original sparse counts |
| Cells, features and design | Every identity, donor, condition and cell-type assignment agrees; 16 samples, eight donors and eight source cell types |
| Source preservation | Every original AnnData slot remains equal; all 28 projected datasets retain stored values/dtypes through annotation |
| QC | Every cell/feature total and nonzero count agrees with the independent sparse reference |
| Pseudobulk | All 124 observed donor/condition/cell-type groups and all 866,103 nonzero aggregate values agree |
| Native reconstruction | Projection, count store and aggregate replay pass |
| Physical Mac mini | Complete count-store and aggregate payload hashes match the laptop; both reconstruct successfully |
| Legacy regression | Twelve positive cases and nine controlled rejections pass under v2 |
| Modern regression | All ten historical unique-selection bundles retain exact original source/plan/output/report bytes |

There are 124 observed groups out of 128 possible combinations; absent
combinations are not fabricated. The independent checker builds group membership
from the original annotations and checks every native membership list before
comparing sparse group sums. It never materializes a dense cells × genes matrix.

On the physical M4 Pro Mac mini, count-store publication/reconstruction took
1.26/0.80 seconds with peak resident memory 59,637,760/62,357,504 bytes.
Aggregation publication/reconstruction took 1.68/1.68 seconds with peaks
207,306,752/221,626,368 bytes. These are recorded native-owner executions on a
shared host, not a controlled speed comparison with scverse or a Metal benchmark.

## Reproduction and provenance

```sh
NUMIVIVO_HDF5_LIBRARY=/path/to/libhdf5.dylib python check_legacy_count_route.py \
  --binary /path/to/h5ad-check \
  --source /path/to/original-kang.h5ad \
  --mapping /path/to/frozen-kang-fit.json \
  --out /new/complete-count-route
```

The mapping argument accepts either the existing fit plan's `mapping` field or
an explicit import mapping. Python prepares plans and independently checks the
results; Swift/native HDF5 performs projection, annotation, import and aggregation.
The scoped harness now exposes the existing count-store public owner through
`count-store` and `verify-count-store`, avoiding a new full product rebuild.
The corresponding product commands remain `singlecell-h5ad-store` and
`singlecell-count-store-verify`.

The tested executable SHA256 is
`20b91b15b9e3c784fd6d441d88152ccb99375a2a75673d996007a73c1fd4b303`.
The scoped receipt implementation tag remains zero; the actual binary, HDF5
library and all 89 compiled owner source files are separately hash-bound.
The original legacy file is retained in the projection bundle; annotation's
journal binds that projected source; count-store and aggregation receipts bind
the exact annotated file. Preserving only the last count file would omit this
provenance chain.

[evidence/2026-09-11-legacy-count-route](evidence/2026-09-11-legacy-count-route)
retains source, plans, complete reports, receipts, independent checks, execution
logs and failures. Large original/projected/annotated H5AD files and binary
counts remain at manifest-listed external paths with exact hashes. Verified
remote duplicates are removable only against those retained copies; the cleanup
manifest records exact restoration mappings.

This qualifies analytical input preservation and arithmetic. It introduces no
new DE fit, held-out prediction, independent biological replication, calibrated
intervals or phenotype result. The broader [prediction assessment](../../../../Documentation/BiologicalPrediction.md)
and [development goal](../../../../Documentation/SingleCellInteroperability.md)
retain their outstanding requirements.
