# Native Visium directory import

`multiassay-visium-import` now consumes an original filtered Visium RNA HDF5
matrix and its Space Ranger position CSV directly. No Python count/coordinate
conversion or intermediate H5MU is required. It uses the existing native 10x
count owner and common multi-assay representation, with observations typed as
spots and a declared full-resolution `[x,y]` pixel frame.

```sh
export NUMIVIVO_HDF5_LIBRARY=/absolute/path/libhdf5.dylib
numivivo multiassay-visium-import /path/to/outs \
  --plan lymph-node-plan.json --output /path/to/new-bundle
numivivo multiassay-verify /path/to/new-bundle
```

The [example plan](lymph-node-plan.json) describes the original complete Human
Lymph Node release. Replace sample, organism, assay namespace, assembly and
source/frame descriptions with the identities of your actual experiment. The
input directory must contain `filtered_feature_bc_matrix.h5` and either
`spatial/tissue_positions_list.csv` (`legacyHeaderlessCSV`) or
`spatial/tissue_positions.csv` (`csvWithHeader`), selected explicitly in the plan.
The latter requires the exact six-column Space Ranger header. This path supports
one RNA UMI assay; it does not silently pick a file when both formats exist.

Coordinates follow the [publisher's definition](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs):
pixel column is x and pixel row is y. The importer validates the six fields,
unique barcodes, binary tissue membership and exact nonnegative integer
coordinates. The filtered count barcode set must equal the complete in-tissue
position set. It joins by barcode and preserves count order; missing, extra,
duplicate or outside-tissue count rows reject. There is no implicit spot or gene
filtering. The position parser is limited to 16 MiB and the existing observation
limit; resident assay limits, including 32 million nonzeros, still apply.

The bundle retains `original.h5`, the exact selected `positions.csv`, the plan,
`dataset.json`, `dataset.h5mu` and a receipt binding both source files. Native
verification reconstructs from these immutable snapshots. The position snapshot
retains all array indices and outside-tissue rows; the typed dataset exposes
selected spots and full-resolution positions. Images and scale-factor files are
not copied into the bundle or interpreted. Image registration, micrometer
calibration, segmentation, deconvolution and tissue prediction are not supplied.
Matrix Market input, Visium HD Parquet/segmented formats, multiple sections and
other assay types require additional explicit support.

## Complete real-data qualification, 2026-09-11

The benchmark reuses the [original Space Ranger 1.0.0 Human Lymph Node release](../SPATIAL.md),
including every filtered spot and gene. Its 4,992 original position rows include
953 outside-tissue rows. A second structural input adds the documented modern
header to those same rows; it is not another biological experiment or a native
newer-release validation.

| Check | Result |
| --- | --- |
| Spots / genes | All 4,039 / 33,538, exact identities and order |
| Nonzero counts / RNA UMIs | All 21,210,753 / 70,070,300, exact |
| Pixel coordinates | Every spot matches source and independent Scanpy |
| Native import and reconstruction | Both CSV formats pass all four product commands |
| Previous qualified dataset and H5MU | Both formats reproduce every byte |
| Swift tests | 13 tests in three suites pass, with native HDF5 enabled |
| Product lifecycle/regressions | 16 commands meet their declared success/rejection outcome |
| Nonspatial 10x compatibility | Previous executable's dataset, H5MU, plan and source bytes preserved |

The lifecycle checks cover deterministic repeat publication, reconstruction after
an original source changes, position tampering, existing destinations, malformed
membership, duplicate barcodes, ambiguous format, wrong units and symlink source
rejection. Rejected imports leave no destination or staging directory. Existing
nonspatial receipts omit the new optional position fingerprint.

On the physical M4 Pro Mac mini (24 GB, macOS 26.6, Swift 6.3.3), import took
8.98/9.18 seconds and verification 9.50/8.77 seconds for legacy/header input.
Peak resident memory was 3,109,781,504 bytes. These are individual operational
measurements, not a matched speed benchmark, out-of-core execution or a GPU claim.
The initial test harness missed the framework and macro plugin search paths;
both failed logs remain beside the successful run. Reference deprecation warnings
also remain. No full-package build or new biological prediction is claimed.

Native executable SHA-256:
`1508dfefb93c695b54e5f2802e2aabbfbf685f005ea95feec2137cff7cc42f7b`.
HDF5 runtime SHA-256:
`a00ffbf8ab94ad81f67231a1ae01df748689e1c35a3615a57f88f4710b8d213e`.
Input freeze:
`f9cf26a66e2ec9c44f8f4e9ad7a17c3d3e539a678e76f9dd13a7dca5804b03e2`.

## Reproduction and retained evidence

From the repository root, build the actual scoped product CLI:

```sh
bash Tools/Omics/H5AD/build.sh /path/to/runtime --with-cli
bash Tools/Omics/Multimodal/Visium/test.sh /path/to/runtime
python Tools/Omics/Multimodal/Visium/prepare.py \
  --prior Tools/Omics/Multimodal/evidence/2026-09-09-spatial --out /path/to/inputs
python Tools/Omics/Multimodal/Visium/run.py --inputs /path/to/inputs \
  --binary /path/to/runtime/numivivo-omics --out /path/to/native
python Tools/Omics/Multimodal/Visium/check.py --inputs /path/to/inputs \
  --native /path/to/native \
  --prior Tools/Omics/Multimodal/evidence/2026-09-09-spatial --out /path/to/checks
```

The independent checker uses h5py 3.16.0, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.5
and Scanpy 1.12.4. Use a separate environment; native import requires no Python.
`prepare_regression.py` creates the declared structural fixture; `regression.py`
compares it with an explicit previous executable. The recipes accept new output
directories and preserve failed runs.

[Compact evidence](evidence/2026-09-11/manifest.json) retains executed sources,
plans, checks, logs and fingerprints. The complete source files remain in the
existing spatial evidence directory. Full native bundles are retained on both
hosts under `numivivo-visium-native-20260911/native/{legacy,header}.tar.gz`.
Their member sizes and SHA-256 values are in the archived execution record.
Together they preserve 1,093,245,305 decoded bytes in 251,660,386 stored bytes.
Every member was checked before removing only this run's finished scratch
copies, with no open handles. This reduces retained size by 841,584,919 bytes;
it is not a claim that total free disk increased relative to the start of work.

To restore a bundle, verify its recorded archive hash, extract its six regular
members into a new directory, verify each recorded member hash and run
`multiassay-verify` with the recorded executable. Receipts bind to that executable;
a different build must produce its own receipt. Spatial analytical methods,
biological ATAC/spatial validation and joint multimodal prediction remain open.
