# Native ragged-array projection and LZF decoding

Native `singlecell-h5ad-project` now preserves AnnData `awkward-array` 0.1.0 fields in `obsm`, `varm` and `raw/varm`. Filtering, reordering, repeated indices and chained projections retain nested lists, records, missing values, unions, integer precision and boolean/string payloads. This extends projection of existing arrays; it does not add a typed API for creating arbitrary ragged annotations or an AIRR analysis model.

The implementation retains the original typed buffers and composes the outer selection. Indexed, option, masked and union roots receive the appropriate index/mask transformation; other roots receive an index wrapper. Raw's independent feature axis stays intact. Unreachable payload remains stored, so this is **not storage compaction or removal of unselected data from the file**. Existing work/output limits apply, with bounded numeric buffer vectors and a 16 KiB structural-form limit. Nested payload is preserved opaquely; this is not a complete validator of every possible malformed Awkward form. Legacy ragged compound datasets remain unsupported.

The storage contract follows [AnnData's Awkward encoding](https://anndata.readthedocs.io/en/stable/fileformat-prose.html#awkwardarrays). Scverse currently marks this encoding experimental. Native execution does not require Python or Awkward; those libraries provide independent qualification.

## Real receptor data

The complete [Scirpy Wu 2020 3k release](https://scirpy.scverse.org/en/latest/generated/scirpy.datasets.wu2020_3k.html) passed: **3,000 cells, 7,544 receptor chains, 20 fields per record**. Identity projection and reversing all cells followed by two repeated cells both preserve every receptor record, field type, original buffer byte and cell-metadata row, and pass native reconstruction. This is the complete published 3k subset, not the full original study.

The input is `https://exampledata.scverse.org/scirpy/wu2020_3k.h5mu`, verified against Scirpy's registry MD5 `12c57c790f8a403751304c9de5a18cbf` and SHA-256 `a28195a12c9758ac738d3583626694020059d8c33691354ccfa6403110ced566`. The checker mechanically copies the complete embedded `mod/airr` AnnData group to a standalone H5AD, preserving its compressed datasets. No cells or receptor fields are removed or synthesized. This does not qualify native H5MU import of AIRR data or the separate GEX modality.

The original LZF-compressed release exposed a missing native decoder. Native HDF5 loading now registers a decoder for filter 32000 when the runtime lacks one. It supports both current size hints and older files without them, grows only on insufficient output capacity, rejects malformed back-references and limits each decoded chunk to 64 MiB. Existing registered filters take precedence. No encoder, Python runtime or external filter plugin is introduced. The [h5py filter](https://github.com/h5py/h5py/blob/master/lzf/lzf_filter.c), [LZF stream decoder](https://github.com/h5py/h5py/blob/master/lzf/lzf/lzf_d.c) and [HDF5 filter ABI](https://github.com/HDFGroup/hdf5/blob/hdf5_1.14.6/src/H5Zdevelop.h) document the format and interface.

## Qualification

The physical Mac mini's scoped Omics library and actual product router built successfully. Final executable SHA-256: `ff922a3911e08e8cf26de64bc9ab1167cb5303fe376764022e4d20ed29c717b8`. AnnData 0.13.3.post0 and Awkward 2.13.0 independently read and compare the native outputs.

- [Ragged conformance](../check_awkward.py): 12 command cases, including every supported root class, empty roots, repeated/chained selections, malformed root buffers and resource-limit rejection.
- [Real-data check](../check_awkward_wu2020.py): all records and original buffers, both projections and reconstruction.
- Existing annotation suite: 26 checks passed. Tensor annotation/projection suite: eight checks passed on the final binary.
- [LZF controls](../LZFControls.swift): AddressSanitizer passes literal, overlapping and long back-references, malformed/truncated streams, capacity bounds, 20,000 deterministic byte streams and 24 independent h5py-compressed chunks. The [generator](../prepare_lzf_controls.py) reproduces those chunks and extracts the unchanged error declaration from the owning source; the actual full `VivoHDF5.swift` is compiled with the sanitizer.

`evidence.tar.gz` retains reports, plans, source hashes, checker versions, logs and the artifact manifest. The initial unsupported-encoding failure, missing-filter failure, no-size-hint decoder failure and first compile errors remain retained. Full source/output H5ADs, upstream H5MU, binaries and compressed test chunks remain under `/Users/n/numivivo-h5ad-ragged-20260912`; they are hash-bound in the manifest rather than duplicated in Git. `qualified-sources.json` includes the filter ABI declaration and module map. The build script now includes those declarations in future source manifests.

This proves the tested data-preservation route. It does not establish receptor specificity, immune recognition, multimodal biological preservation or clinical efficacy.
