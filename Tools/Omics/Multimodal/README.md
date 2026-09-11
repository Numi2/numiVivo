# Native multi-assay count data and CITE-seq interoperability

The common model now represents separate feature spaces over shared observation
identities. `VivoMultiAssayDataset` owns sample/donor metadata, observations and
spatial frames; each `VivoAssaySpace` owns its namespace, assay kind, count unit,
feature definitions, local-to-global observation map and exact sparse counts.
RNA, antibody capture, accessibility and guide capture are distinct assay kinds.
An assay can cover a subset of observations in an independent row order.
A missing row means unmeasured; a measured sparse zero remains zero. Duplicate
feature IDs are allowed across spaces, but not within a space.

This is a count-assay foundation, not a joint latent model or a complete
multimodal analysis suite. Continuous protein/metabolite measurements, joint
RNA/ATAC analysis, WNN integration and biological spatial/ATAC qualification remain open. Current native single-cell algorithms retain their
own RNA/count interfaces; this change does not silently apply RNA normalization
or differential expression to antibody or accessibility features.

## Product commands

```sh
numivivo multiassay-10x-import original.h5 --plan pbmc5k-plan.json --output bundle
numivivo multiassay-h5mu-import source.h5mu --plan pbmc5k-h5mu-plan.json --output imported
numivivo multiassay-verify bundle
numivivo multiassay-h5mu-write bundle/dataset.json --output exported.h5mu
```

The import retains the exact `original.h5`, explicit plan, canonical typed
`dataset.json`, native `dataset.h5mu` and implementation-bound receipt. Verification
reconstructs the count model from retained source and plan, then regenerates H5MU
and compares bytes. Rehashing an edited dataset does not bypass reconstruction.
New destinations are required; failed imports remove their private staging data.

The [10x feature-barcode HDF5 format](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices)
stores features-by-barcodes in CSC. The native reader consumes one barcode's
entries using slices of at most 65,536 values, merges duplicate coordinates,
discards explicit zeros, and builds canonical per-assay CSR matrices. Every
source feature type must have an explicit compatible mapping; no assay is silently
discarded. Nonempty source genome annotations must match the declared assembly.
Counts remain UInt64; no float conversion occurs. Type mapping, namespace, units,
sample identity and evidence are supplied explicitly in the plan.

H5MU export follows the [MuData layout and one-based observation/feature maps](https://mudata.readthedocs.io/stable/io/spec.html).
Each modality is an AnnData containing its own matrix and feature IDs. Global
observation IDs are unique generated indices; exact source sample/barcode pairs
remain separate columns. Observation maps preserve missing assays and row order.
Complete typed metadata, including spatial frames and genomic intervals, is
retained in `uns['numivivo']` JSON. Spatial coordinates are also exported as numeric
`obsm` arrays with explicit frame-to-path metadata; see [spatial interchange](SPATIAL.md). MuData emits a
creator warning for this native HDF5 producer; independent reads and roundtrips
are explicitly checked. The original 10x file retains any source-specific tags
outside the typed count projection.

Accessibility features require an assembly and explicit zero-based half-open
intervals. The 10x Peaks mapping parses `contig:start-end`; the caller declares
fragment, read or cut-site count units, according to the source measurement. Spatial observations can be cells or spots, with
finite 2D/3D positions in named frames declaring axes and pixel/micrometer units.
These representations are structurally tested, not biological qualifications.

## Native H5MU count import

`multiassay-h5mu-import` supports MuData 0.1.0 with shared observations (axis 0),
AnnData 0.1.0 modalities and dataframe 0.2.0 axes. Each modality must be mapped
explicitly to a feature namespace, assay kind, count unit and optional assembly.
Select `X` or an explicit `layers/<name>` count matrix independently per modality.
CSR, CSC and dense numeric layouts use the same exact count decoder as H5AD;
normalized fractional values reject if selected. `raw/X` requires a different
feature-axis model and is not accepted by this importer.

Global sample/barcode/group/kind columns define observations. Selected local
columns must agree wherever present. One-based `obsmap` and `varmap` values are
validated against global/local names: zero means absence, each local row/feature
must occur once, and feature spaces cannot overlap in the global variable map.
Explicit source namespaces, assay kinds and genome annotations cannot contradict
the plan. String, categorical and nullable-string identity columns are supported.

An optional spatial mapping selects a numeric `obsm/<name>` array and declares
its frame, axes and units. Entirely NaN rows mean missing positions; partially
missing or infinite positions reject. Native export retains these positions in
typed `uns` metadata and standard numeric arrays as described above. Accessibility mappings require the explicit convention
`contig:start-end:zero-based-half-open`, an assembly, and fragment/read/cut-site units.

The import is a declared count/identity projection. Unselected annotations and
other analysis objects remain in the exact retained source file. Every source
format is retained as `original.h5`; receipts select the appropriate reconstruction
reader. Source and projection hashes remain distinct.

## Public real-data input

Use the complete filtered
[10x PBMC5k TotalSeq-B dataset](https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-with-cell-surface-proteins-v-3-chemistry-3-1-standard-3-1-0),
Cell Ranger 3.1.0, published 2019-07-24, licensed CC BY 4.0 by 10x Genomics.
It is one healthy human donor; it is not a multi-donor or replicated perturbation
experiment. The page describes a 31-antibody panel; the actual full source
feature table contains 32 antibody entries, all retained without interpreting
or discarding that discrepancy.

```sh
curl -fL 'https://cf.10xgenomics.com/samples/cell-exp/3.1.0/5k_pbmc_protein_v3/5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5' -o original.h5
python check_citeseq.py --binary /absolute/path/numivivo --source original.h5 --out qualification
```

Pinned bytes: 17,113,138. SHA-256:
`3b290ad9605b96974c9c16e5ae3427e5e5c496a66e55223df135b388b6d61417`.
The full source contains 5,247 cells, 33,538 RNA features and 32 antibody features;
RNA has 9,250,566 nonzeros and antibody capture has 155,429. No cells or features
are filtered for the comparison. The checker uses h5py/SciPy as the exact-count
oracle, Scanpy's full feature reader, and MuData for native H5MU read/roundtrip.

## Execution and bounds

The current model is resident after import: at most 100,000 global observations,
16 assays, 200,000 total features and 32 million total retained/source entries.
One source barcode's entries, indices and metadata are also resident. H5MU
imports additionally hold transient row dictionaries before final CSR encoding;
dense scan work is limited to 500 million elements across all modalities. Sources
and H5MU files are capped at 1 GiB, canonical dataset JSON at 512 MiB. This is not
million-cell, memory-mapped multimodal execution or a CPU/GPU speed claim.

HDF5 remains an optional dynamic runtime dependency. Pure schema tests need no
HDF5. The native HDF5 structural test is opt-in through `NUMIVIVO_TEST_HDF5=1`;
set `NUMIVIVO_HDF5_LIBRARY` when using a nonstandard native library location.
The qualification used an isolated Mac mini HDF5 runtime copied from the laptop's
HDF5 2.2.0 and libaec libraries, with local loader paths and ad-hoc signatures.
No global package manager or runtime was installed on the Mac mini.

## Qualification on 2026-09-09

The release build passed and all 58 Swift tests in 16 suites passed, including
native HDF5 with the opt-in flag enabled. Two complete product imports produce
byte-identical datasets, H5MU exports and receipts. The scoped and full product
builds also produce identical dataset/H5MU bytes. The 23 product CLI cases include
17 expected rejections; six of those rejection cases and two positive commands
are explicitly synthetic RNA/ATAC structural controls. They supplement the
complete real CITE-seq comparison rather than establishing biological ATAC quality.

All 32,173,181 RNA UMIs and 16,081,143 antibody UMIs, all nonzero coordinates,
cell/feature axes, per-cell totals and assay memberships agree with the original
h5py/SciPy matrix and Scanpy's full feature reader. MuData 0.4.1 reads both exported
modalities and preserves their counts and alignment through its own roundtrip.
Missing-assay row maps, spatial frame metadata and UInt64 values above 2^53 are
verified separately with structural controls. Malformed row maps, frames, units,
feature mappings, assemblies, counts, peak intervals and sparse indices reject;
source reconstruction also rejects deliberately rehashed count changes.

The initial compile rejected two missing `try` expressions. After repair, an
initial test run failed because the Mac mini lacked HDF5. Both failures are
retained in evidence; the final run explicitly enabled HDF5 and passed all tests.
Product/source/test fingerprints and native library hashes are retained in
`evidence/2026-09-09/qualification.json`. The HDF5 and library-only dependencies
are not bundled into the NumiVivo product by this change.

To reconstruct the archived verified bundle, copy `original.h5`, `plan.json`
and `receipt.json` from the evidence directory, then decompress `dataset.json.gz`
and `dataset.h5mu.gz` into that new bundle. Receipts bind to the recorded executable;
a rebuilt executable can create and verify its own fresh bundle. The original
source is small enough to retain directly here under the publisher's CC BY 4.0
license, with the attribution and source link above.

## H5MU import qualification on 2026-09-09

All 61 selected Swift tests in 17 suites passed with native HDF5 enabled; the
release and scoped builds passed. The full 5,247-cell CITE-seq dataset imported
from native CSR H5MU, independent MuData CSR H5MU, and mixed RNA CSC / protein
dense count layers produces exactly the previously qualified dataset bytes.
The mixed-layout source deliberately contains fractional values in unselected X.
Repeated imports produce identical source, plan, dataset, H5MU and receipt bytes.

`check_h5mu.py` passed 32 product commands, including 23 expected rejections.
Structural controls preserve reordered partial assay rows, exact UInt64 values
above 2^53, spots and explicit spatial positions. Invalid maps, metadata conflicts,
count values, paths, spatial positions and rehashed output tampering reject.
`check_h5mu_atac.py` adds six product commands for synthetic peak mapping,
reconstruction and four expected convention/assembly/unit rejections.

Shared-reader regressions passed the AnnData interoperability suite, 25 streaming
checks, complete real PBMC3K preservation/count/QC/normalization comparisons, and
the existing full CITE-seq 10x checker. These are interoperability and numerical
checks, not qualification of joint multimodal inference or real spatial/ATAC biology.

```sh
python check_h5mu.py --binary /absolute/path/numivivo --citeseq qualification --out h5mu-qualification
python check_h5mu_atac.py --binary /absolute/path/numivivo --out h5mu-atac-controls
```

Evidence is in `evidence/2026-09-09-h5mu`. The first Python checker stopped when
its mutation helper assumed categorical storage for a nullable-string column.
That failed log is retained. After correcting the helper, the complete checker
passed; no native implementation change was needed. MuData's native-producer and
cross-space duplicate-variable-name warnings remain in the logs.

The [complete paired RNA/ATAC benchmark](MULTIOME.md) documents the separate
Cell Ranger ARC cut-site unit and the full-source qualification protocol.

## Native Visium directory input

The [Visium importer](Visium/README.md) now accepts original filtered RNA HDF5
and explicit legacy/header Space Ranger position CSVs. It preserves full-resolution
pixel positions by barcode and binds both original files in its native receipt.
The full 4,039-spot source passes exact counts/coordinates and reconstruction,
without the earlier external H5MU preparation. This is resident interchange;
HD formats, image processing and spatial biological prediction remain open.
