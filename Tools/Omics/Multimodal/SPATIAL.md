# Spatial count and coordinate interchange

Native MuData export now writes each authored spatial frame as a numeric global
`obsm` array. The first available name is `spatial`, followed by `spatial_1`,
`spatial_2`, and so on. Names already used by modality membership masks are skipped.
`uns['numivivo']` records the exact frame-ID-to-array-path mapping as `spatialArrays`,
alongside frame axes, pixel/micrometer units and source descriptions. Frame IDs
are metadata, never filesystem/HDF5 paths. Observations in a different frame or
without a position have entirely NaN rows in that frame's array.

A single-frame export with the usual `obsm/spatial` name can be directly reimported
with the same explicit spatial mapping. Multiple frames are exported separately;
the H5MU count importer currently selects one spatial array per plan. It does not
merge coordinate frames or infer physical transformations. Sources without spatial
frames retain their previous native H5MU representation.

## Public input and coordinate convention

The benchmark uses the complete filtered [10x Human Lymph Node Visium dataset,
Space Ranger 1.0.0](https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-0-0),
licensed CC BY 4.0. It includes 4,039 tissue spots and 33,538 genes from one section.
No spots or genes are filtered by the benchmark. A spot is represented as a spot,
not an inferred single cell. The later Space Ranger 1.1 release is a different
source and is not substituted here.

The [Space Ranger spatial-output definition](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs)
identifies `pxl_col_in_fullres` as x and `pxl_row_in_fullres` as y. The plan declares
`[x,y]` axes in full-resolution **pixels**. The checker joins these coordinates by
exact barcode, verifies that the retained barcode set equals the in-tissue set,
and compares the coordinates independently with Scanpy's Visium reader. There
is no conversion to micrometers or calibration inferred from nominal spot size.

## Reproduce

```sh
curl -fL 'https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5' -o original.h5
curl -fL 'https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_spatial.tar.gz' -o spatial.tar.gz
python check_spatial.py --binary /absolute/path/numivivo --source original.h5 --spatial spatial.tar.gz --out qualification
```

| Input | Bytes | SHA-256 |
| --- | ---: | --- |
| Count HDF5 | 27,416,993 | `eb883e48d9aa935d7959600153ad7ba8c7f1748302fcf56eef4c6999574bb37d` |
| Spatial archive | 8,733,712 | `064f7e62e43a705730c30911911940dbbe583ebc6b743dd2dc52eba5037a2ef5` |

Reference preparation constructs MuData from the original count matrix and
barcode-joined positions. This earlier Python preparation is explicit. A subsequent
[native Visium directory importer](Visium/README.md) now reads the original
filtered RNA HDF5 and selected position CSV directly; this section retains the
original H5MU-based experiment. The native path imports H5MU, preserves
spots/counts/coordinates, exports standard arrays and verifies source reconstruction.
The exact original count and spatial files are retained in the evidence archive.

The checker compares every count and stable feature identity with h5py/SciPy and
Scanpy, checks all coordinates, reimports the native export and requires identical
native dataset/projection bytes. A standalone image overlay uses the supplied
low-resolution scale factor solely for display. A separate synthetic control
checks 2D/3D frames, missing rows and collisions with an assay named `spatial`.

This is interchange evidence. It does not qualify cell segmentation, image
registration accuracy, spatially variable genes, spot deconvolution, physical
length calibration or tissue prediction. The count model remains resident.

## Qualified result on 2026-09-09

The release build and all 63 selected Swift tests in 17 suites passed with native
HDF5 enabled. Six spatial product commands passed, covering real-data import,
source reconstruction, repeatability, native reimport, reimport reconstruction
and a synthetic multiple-frame/name-collision control. A separate full Multiome
re-export confirmed that nonspatial H5MU bytes remain unchanged.

All 21,210,753 nonzero entries, 70,070,300 RNA UMIs, 4,039 spot coordinates and
33,538 stable gene identities agree with the source and independent readers.
Repeated native bundles and native spatial reimport are byte-identical. The
[source/native coordinate overlay](evidence/2026-09-09-spatial/coordinate-overlay.png)
was visually inspected; it demonstrates matched plotting coordinates, not an
independent assessment of the publisher's image registration.

Maximum resident memory across the spatial commands was 4,629,364,736 bytes
(4.63 GB) on Apple M4/macOS 26.6. The release executable was built on the Mac mini
and used an isolated dynamic HDF5 2.2.0 runtime. This is CPU execution, without
Metal acceleration or an out-of-core claim.

Evidence is in `evidence/2026-09-09-spatial`, including original source HDF5 and
spatial archive, prepared H5MU snapshot, native dataset/projection, plans, receipts,
logs, runtime/source fingerprints and overlay. To restore the verified native
bundle, decompress its `original.h5.gz`, `dataset.json.gz` and `dataset.h5mu.gz`,
and copy its plan and receipt. Receipts bind to the recorded executable. A new
build must generate its own receipt.

The pinned Scanpy version warns that `read_visium` is deprecated; this reference
check still passed. Duplicate gene display-name and MuData creator warnings are
retained. Stable feature IDs were not renamed. The subsequent [native directory qualification](Visium/README.md) covers
filtered RNA HDF5 plus legacy/header CSV. Multiple-frame import in one plan,
Visium HD formats and spatial analytical methods remain open.
