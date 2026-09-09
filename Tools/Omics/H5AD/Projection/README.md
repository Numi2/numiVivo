# Native AnnData axis projection

`singlecell-h5ad-project` creates a new, source-bound H5AD bundle after selecting
or reordering unique cell and feature indices. It transforms every supported
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
that complete axis; an empty array creates an empty axis. Duplicate indices and
out-of-range indices are rejected. Source identity is checked against the copied
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
| `uns` | Copied without inferred alignment | Copied without inferred alignment |

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
axis. Sparse structural indices/offsets are written as int64 to represent remapped
positions without narrowing. Compression, chunk layout and hard-link alias
identity are not promised; values behind aliases remain correct in their
respective aligned or unstructured contexts.

Unknown aligned encodings, unknown root fields, external/soft links, external or
virtual dataset storage, object references, null dataspaces requiring projection,
first-axis/matrix shape mismatches, and duplicate selection indices fail explicitly. General legacy
unversioned formats, ragged/Awkward arrays, structured record encodings and
axis duplication remain outside this operation. No unsupported slot is silently
removed to produce a successful result.

## Storage and replay

The bundle contains `original.h5ad`, `projected.h5ad`, `plan.json`, `report.json`
and `receipt.json`. Source, plan, result and report hashes plus the implementation
fingerprint are bound by the receipt. Replay snapshots the original, repeats the
projection and checks both the report and exact H5AD output bytes. New HDF5 objects
have timestamp tracking disabled so replay does not depend on the wall clock.
Failures remove the staging directory; existing destinations are not overwritten.

Source and output files are bounded to 1 GiB. Dataset allocations are checked
before creation, variable storage is reserved during block transfers, and file
size is checked during writes/copies. Axes are bounded to two million entries;
each sparse input has at most 100 million stored entries. The default work limit
is 500 million element visits, configurable up to two billion via
`maximumElementVisits`. Dense output visits and sparse index/data passes share
this allowance. Empty arrays consume zero element visits.

Transfers gather at most 65,536 values at once with bounded numeric and variable
buffers. Dense source arrays are read in blocks, and sparse inputs remain sparse;
there is no dense cells × genes intermediate. Axis maps, sparse offsets and
per-axis counts remain resident. These bounds and tests do not establish
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
empty selections. Rejections cover wrong source hashes, duplicate/out-of-range
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
