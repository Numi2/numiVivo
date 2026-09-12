# Native AnnData axis projection

`singlecell-h5ad-project` creates a new, source-bound H5AD bundle after selecting
or reordering cell and feature indices, including repeated selections. It transforms every supported
aligned slot, preserving stored values and datatypes instead of reconstructing
only an expression matrix. The original file remains in the bundle.

```sh
numivivo singlecell-h5ad-project source.h5ad --plan projection.json --output projected
numivivo singlecell-h5ad-project-verify projected
```

For example, generate a plan from the actual source bytes:

```python
import hashlib, json
with open("source.h5ad", "rb") as source:
    fingerprint = list(hashlib.file_digest(source, "sha256").digest())
with open("projection.json", "w") as out:
    json.dump({
        "schemaVersion": 1,
        "source": {"bytes": fingerprint},
        "provenance": "Explicit analysis cohort and feature selection",
        "observationIndices": [5, 1, 0],
        "featureIndices": [8, 2, 3]
    }, out)
```

Indices refer to the original zero-based axis order. Omitted/null indices keep
that complete axis; an empty array creates an empty axis. Repeated indices repeat
the corresponding aligned values in the requested order. Negative and out-of-range
indices are rejected. Source identity is checked against the copied
snapshot before projection. Python above only prepares JSON; execution and replay
use Swift and native HDF5 without Python/scverse.

## Alignment rules

| Slot | Cell selection | Feature selection |
| --- | --- | --- |
| `X`, each `layers` matrix | Rows | Columns |
| `obs` dataframe | Rows, index and every column | None |
| `var` dataframe | None | Rows, index and every column |
| `obsm` arrays/sparse arrays/dataframes | First axis | None |
| `varm` arrays/sparse arrays/dataframes | None | First axis |
| `obsp` arrays/sparse arrays | First two axes | None |
| `varp` arrays/sparse arrays | None | First two axes |
| `raw/X` | Rows | Retains the independent raw feature axis |
| `raw/var`, `raw/varm` | None | Retained in full |
| `uns` | Copied without inferred alignment; eligible legacy categories move into obs/var | Same |

These rules follow the [AnnData storage specification](https://anndata.readthedocs.io/en/stable/fileformat-prose.html)
and [raw slicing semantics](https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.raw.html).
Existing analysis metadata in `uns` is preserved; projection does not recompute
neighbors, normalize induced subgraphs, refit PCA or regenerate differential
expression. The operation records index selection and source provenance rather
than claiming those analyses remain valid for a changed cohort.

Supported aligned encodings are dense numeric/string arrays, CSR/CSC matrices,
dataframes, categoricals, nullable integer/boolean/string columns, and raw groups
under their current versioned AnnData encodings. Arrays can have up to eight
axes; matrix/layer inputs must have two. Unused categories and their order are
retained, including during empty selection. The independent AnnData comparison
therefore sets `remove_unused_categories=False` explicitly.

Numeric data never pass through Double. This preserves integer values above
2^53, endian-specific floats, NaN/Inf payloads, boolean enums, and complex numeric
arrays where encoded as ordinary arrays. Sparse values retain their dtype,
explicit zeros, duplicate entries and source order within each selected major
axis. When a minor index is selected repeatedly, each stored source entry expands
to its requested destination positions in ascending output order. Duplicate
source entries remain distinct; repeated major axes repeat their whole stored
sequence. Sparse structural indices/offsets are written as int64 to represent remapped
positions without narrowing. Compression, chunk layout and hard-link alias
identity are not promised; values behind aliases remain correct in their
respective aligned or unstructured contexts.

Unknown aligned encodings, unknown root fields, external/soft links, external or
virtual dataset storage, object references, null dataspaces requiring projection,
first-axis/matrix shape mismatches fail explicitly. The legacy route below
supports specified compound dataframe and embedding representations. Other
legacy formats, ragged/Awkward arrays and arbitrary nested record encodings
remain outside this operation. No unsupported slot is silently removed.

## Storage and replay

The bundle contains `original.h5ad`, `projected.h5ad`, `plan.json`, `report.json`
and `receipt.json`. Source, plan, result and report hashes plus the implementation
fingerprint are bound by the receipt. Replay snapshots the original, repeats the
projection and checks both the report and exact H5AD output bytes. New HDF5 objects
have timestamp tracking disabled so replay does not depend on the wall clock.
Failures remove the staging directory; existing destinations are not overwritten.

Source snapshots use the existing streamed owner's 64 GiB input bound. This
experiment does not qualify a 64 GiB source. Output defaults to 1 GiB; an explicit
`maximumOutputBytes` plan field permits a positive allowance up to 8 GiB.
Omitting that field preserves historical plan bytes and the default bound.
Dataset allocations are checked before creation, variable storage is reserved
during block transfers, and actual file size is checked during writes/copies.
The allowance belongs to each projection's HDF5 handle, not global mutable state.
Axes remain bounded to two million entries. Sparse input size is admitted against
the remaining work allowance before integer conversion and scanning, replacing
the former fixed 100-million-entry cap. The default work limit is 500 million
element visits, configurable up to two billion via `maximumElementVisits`.
Dense output visits and sparse index/data passes share this allowance. Empty
arrays consume zero element visits.

The [complete Replogle UPR assay partition](../../PerturbationPrediction/Replogle2020/README.md)
qualifies this extension on 129,839,577 original entries. The RNA output explicitly
requests 2 GiB; the guide output uses the unchanged default. Both reconstruct
exactly, while the full RNA request with default storage or insufficient work
rejects without publication. This is input/assay qualification, not a prediction
result.

Transfers gather at most 65,536 values at once with bounded numeric and variable
buffers. Dense source arrays are read in blocks, and sparse inputs remain sparse;
there is no dense cells × genes intermediate. Axis maps, sparse offsets and
per-axis counts remain resident. Repeated sparse expansion is counted against
the shared work allowance before output allocation, with checked accumulation
when both axes repeat. Expanded transfers still hold at most 65,536 values.
These bounds and tests do not establish
million-cell scalability or a throughput advantage over scverse.

## Independent verification

`check.py` exercises storage fixtures and explicitly supplied real H5AD files:

```sh
python check.py --binary /path/to/numivivo --out /new/check-directory \
  --real /path/to/kang/original.h5ad --real /path/to/baron/original.h5ad \
  --real-mapping /path/to/kang/plan.json --real-mapping /path/to/baron/plan.json
```

A mapping argument may contain the original import mapping or a streaming plan
with a `mapping` field. `--scoped` runs the same tests using `H5AD/build.sh`'s
`h5ad-check` executable. This is an independent Python reference environment;
AnnData, h5py, NumPy, pandas and SciPy are not native runtime dependencies.

The reference compares all AnnData slots, raw axes, dataframe dtypes/categories,
stored HDF5 dtypes/attributes, masked payloads, dense bytes and the exact sparse
entry sequence after remapping. It also reimports each real projected file with
the native count importer, checking counts, feature/cell/sample axes and byte-exact
retention of the projected original H5AD.

Real source scopes are all 2,651 Kang B cells × 15,706 genes and all 8,569 Baron
cells × 20,125 genes. Deterministic reversed/subsampled axes yield 1,326 × 5,236 and
4,285 × 6,709 outputs respectively. These are format/interoperability results,
not new biological benchmarks or validated filtered-cohort inference.

Storage fixtures cover CSR, CSC, dense big-endian arrays, nullable columns,
unused/ordered categories, raw with a different feature axis, tensor embeddings,
dataframe embeddings, pairwise matrices, complex values, large unsigned integers,
nonfinite floats, sparse duplicates, hard-link aliases in `uns`, identity and
empty selections. Rejections cover wrong source hashes, negative/out-of-range
indices, work limits, unsupported encodings, misaligned arrays, output allocation
limits, no-overwrite and tampered output.

The evidence archive retains unsuccessful reference assertions (copied `uns`
index dtype and native cell-field naming), the zero-element budget failure, and
corrected runs. Format fixtures establish software behavior only; real-data
interop and the broader scientific roadmap remain distinct.

## Qualified product snapshot, 2026-09-09

The release CLI passed 10 positive projection cases and 10 rejection cases.
Both real projected files reimported into NumiVivo with exact counts, axes and
retained originals: 258,746 nonzeros for Kang and 2,660,360 for Baron. All 44
existing single-cell tests, the scoped H5AD build, the full release build and
transferred-binary signature/hash checks passed. A later replay of the first
format bundle also passed, beyond the initial creation time.

The qualified release SHA256 is
`766c5f05e5c94da60489f1b19025db8c7122efcce48ee300929d704805cea9e2`.
[evidence/2026-09-09](evidence/2026-09-09) includes source hashes, build/test logs,
scoped and product summaries, and compressed original/projected H5AD bundles
with plans, reports and receipts. Decompress the H5AD files into a new directory
before replay. The two real source snapshots are included so replay does not
require recreating their exact HDF5 encoding from a public dataset download.
These files retain their source-study metadata; the test does not assign new
cell labels or establish scientific validity of the selected cohorts.


## Repeated-index qualification, 2026-09-11

Repeated observation/feature indices now work across all supported aligned slots,
including CSR/CSC, dense arrays, nullable/categorical columns, tensor/dataframe
embeddings, pairwise matrices and raw's independent feature axis. Names and
annotations repeat exactly; the projection does not invent new biological cells,
rename identifiers or treat repeated rows as independent experimental replicates.
The analytical count model still requires unique cell identities. A matched
single-cell control imports, while its repeated-cell counterpart is rejected.

The current native owner passes **19 positive and 11 rejection cases**, including
three transfer-boundary cases with more than 65,536 repeated positions and a
sparse expansion beyond the work budget. All stored datatypes, sparse entries,
AnnData slots and replay outputs agree with the independent reference. Seven
repeated-axis bundles also produce identical bytes on the physical Mac mini.
All ten previously qualified unique-selection bundles reproduce exact original
source, output, plan and report bytes, including both earlier real-data cases.

The new real cases retain every original cell and feature, reverse their order,
and append three repeated positions on each axis:

| Source | Complete source shape | Repeated output shape |
| --- | --- | --- |
| Prepared full Kang | 24,673 × 15,706 | 24,676 × 15,709 |
| Full Baron | 8,569 × 20,125 | 8,572 × 20,128 |

The prepared Kang snapshot is the earlier count/identity-verified conversion with
its declared selected metadata. This result preserves every slot of that prepared
file; it does not qualify the legacy downloaded compound-dataframe encoding or
restore annotations omitted by the earlier preparation. The first new test
attempt encountered that unsupported legacy source before real projection; its
completed fixture evidence and failure remain retained. A separate native call
also rejects the legacy source without publishing partial output.

The ordinary unique real projections also reimport into the count model with
exact counts/axes and retained H5AD bytes: 2,494,517 Kang and 2,660,360 Baron
nonzeros. The complete repeated Kang file exceeds the resident count importer's
file bound; file projection's larger allowance does not remove that separate
limit. A small repeated-identity control exercises the identity rejection
independently. The initial checker expected the word `unique`, while the native
owner reported its existing `cell identity, sample reference or group` diagnostic;
the matched control and explicit diagnostic checks now pass.

The Omics/artifact scoped build passed with all 88 compiled source hashes checked.
The actual tested binary SHA-256 is
`cf5e264d0c4e31cf9f6446cd9b391212265fd785a6d10b2e4c3364f292a84be3`.
Its transfer signature and both-host binary identity were verified. Receipt
implementation tags from the scoped driver are retained separately from that
binary identity. The full unrelated release product was not rebuilt in this turn;
these are native public-API interoperability results, not new biological,
throughput or million-cell qualification.


The [repeated-axis evidence archive](evidence/2026-09-11-repeated) contains 465
members and 2,387,932 stored bytes; manifest SHA-256
`864130f2d3d60bf2ae9d3c81576194e296e07996794ba993aa62dbe9272d14e4`.
Format fixtures and their complete native bundles, every real plan/report/receipt,
source snapshots and all check logs are stored. Twenty-nine large real H5AD/count
payloads remain externally retained with exact paths and hashes; no full-real
matrix was sampled for comparison. Restore the manifest's external files when
replaying the complete archived study on another host.

## Original legacy Kang, 2026-09-11

The native projection owner now admits an unversioned AnnData root, legacy
`h5sparse_format`/`h5sparse_shape` CSR/CSC matrices, compound obs/var datasets,
fixed-array compound obsm/varm fields, and dotted `raw.X`/`raw.var`/`raw.varm`.
It produces versioned aligned slots while retaining the original source file.
The first compound dataframe field becomes its index, and all remaining fields
retain their order. Scalar integer, floating, enum and string fields preserve
their stored representations. String columns use nullable string encoding with
an all-false mask so empty projections retain pandas string semantics.

Legacy string category definitions in `uns/*_categories` follow the installed
AnnData reader's migration rules for obs/var. Every original code is checked
before selection; a code at or above the category count leaves the numeric
column and definition in place. Valid codes retain missing `-1`, unused labels
and label order, with `ordered=false`. Successfully migrated definitions leave
`uns`; unrelated unstructured data remains copied. Raw feature annotations do
not inherit obs/var category migration. Contradictory metadata and invalid codes
fail. Nonstring indices, scalar embedding members and
arbitrary nested/ragged records remain unsupported rather than inferred.

`check_legacy.py` compared the complete **original** Kang file, SHA256
`e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830`,
with AnnData 0.13.3.post0. Both the complete **24,673 × 15,706** output and
reversed full axes plus two repeats (**24,675 × 15,708**) pass. Every AnnData
slot, all 17 compound fields, category semantics, numeric bytes and exact sparse
stored-entry sequences agree. Both outputs reconstruct exactly. The original
matrix has **14,184,532 stored entries**; no cells or genes were removed to make
this qualification succeed. This supersedes the earlier prepared-file evidence
only for the specific legacy input-preservation gap, not its historical analyses.

The native legacy suite passes **12 positive cases and nine controlled
rejections**, covering CSR/CSC/dense, full/repeated/empty axes, fixed-array
tensors, dotted raw, shared obs/var categories, scalar string categories and
unsupported/oversized inputs. All **15 current modern-format cases and 11
rejections** pass. Replaying **10 historical unique** and **19 prior repeated-suite**
bundles preserves the original source, plan, output and report bytes exactly.

All twelve legacy cases were also published and reconstructed on the physical
Mac mini, matching the laptop's four payload hashes exactly, including both
complete Kang projections. After that check, 29 newly created remote H5AD
duplicates were removed only after every retained laptop copy and open-handle
check passed. This removed 476,446,015 payload bytes and measured 399,851,520
bytes of additional APFS availability; logs and restoration mappings remain.

A size query in HDF5 2.2.0 crashed on scalar variable strings during development.
The scalar path now uses a custom allocation bound of 16,385 bytes, allowing at
most 16,384 UTF-8 bytes plus the terminator, without that query. Its positive and
oversized-input rejection cases pass. The crash trace, first compile errors,
pandas big-endian slicing limitation, empty-string inference mismatch and
corrected checks remain in the evidence rather than being counted as passes.

Reproduce using the scoped native owner:

```sh
NUMIVIVO_HDF5_LIBRARY=/path/to/libhdf5.dylib python check_legacy.py \
  --binary /path/to/h5ad-check --out /new/legacy-checks \
  --kang /path/to/original-kang.h5ad
```

The tested binary SHA256 is
`2ad4a4872e7da646ebcdbcbca746dfc7273b7b64bf87c29a7aace1dceaab6fef`,
built on the physical M4 Pro Mac mini from 89 hash-verified owner files. This
scoped harness keeps its zero implementation tag in receipts; the executable
and actual compiled sources are bound separately. No full product rebuild or
new biological prediction result is claimed for this change.

[evidence/2026-09-11-legacy](evidence/2026-09-11-legacy) retains source, reference
reader, checks, failures, format artifacts and every real plan/report/receipt.
Large real H5AD files remain externally retained with restoration paths and
SHA256 bindings in the manifest. AnnData is used only by the independent checker;
native projection and reconstruction require no Python/scverse runtime.

## Native analytical handoff, 2026-09-11

The [complete original Kang analytical route](LEGACY_COUNT_ROUTE.md) now passes
projection → source-bound annotation → streaming count store → pseudobulk.
Legacy projection v2 adds a missing dictionary tag to legacy `uns`, preserving
its children so native annotation can use it. The v1 artifact snapshots above
remain valid historical AnnData-preservation evidence; their native annotation
handoff failure is separately retained. Modern projections keep their exact
historical bytes. The linked report gives every count, group, memory measurement,
source identity and remaining biological limit.


## Numeric legacy categories, 2026-09-12

Legacy numeric category labels now migrate with their original datatype and
bytes. Signed and unsigned integers up to 64 bits are checked without floating
conversion; Float32/Float64 labels retain signed zero and infinities. Duplicate
labels (including opposite-signed zero) and NaN labels are rejected before
publication. Missing category codes retain their original meaning.

The isolated local native build passes 24 reference cases over full, repeated
and empty selections. Cases cover Int64 extrema, UInt64 values above 2^63,
big-endian integers, floats, duplicates and nulls. Pandas cannot read big-endian
category indices directly: that case uses byte-order-normalized disposable
reference copies, with a separate exact datatype/byte check on native output.
Four historical string-category cases preserve all four payload hashes exactly
and pass native reconstruction. These are format tests, not biological evidence.

The [retained evidence](evidence/2026-09-12-numeric-categories/manifest.json)
binds the executable, compiled sources, checker, original failure and all outputs.
Run `check_numeric_categories.py --binary H5AD_CHECK --source integer.h5ad
--out NEW_DIRECTORY` with the retained seed file after extracting the archive.
The broader legacy checker also now expects numeric-category migration to succeed.
