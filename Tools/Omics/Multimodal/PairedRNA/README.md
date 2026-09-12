# Paired RNA PCA qualification

The complete RNA assay from the same 10x paired benchmark used for ATAC LSI
now passes native sparse RNA reduction against executed Scanpy and SciPy.
All 2,711 nuclei, 36,601 genes, 5,218,473 nonzero records and 11,786,194 RNA UMIs
are retained as input. RNA identities match every ATAC observation exactly.

The frozen protocol uses library-size normalization to 10,000 over all genes,
log1p, 2,000 variable genes in 20 mean bins, 30 centered PCA components,
a 256-vector basis, seed 7 and native residual tolerance 1e-6. No mitochondrial,
cell-type or outcome selection is applied. The native positive-count/feature
filter retains every nucleus. Mitochondrial annotations are not supplied or used.

The existing native 10x reader, processing and reduction implementation execute
through a scoped harness linked to the actual library built for commit
74a76d141d4f36c0379c7abc470c77a23ac0ba6b. The harness uses test visibility for
the existing reduction owner; no alternate PCA implementation is substituted.
The explicit resident sparse limit is six million nonzero RNA records; this
does not qualify million-cell execution or a new product CLI.

| Comparison | Result |
| --- | ---: |
| Scanpy 1.12.4 Seurat-flavor HVG identities | All 2,000 exact |
| Relative singular-value error | 2.34e-15 maximum |
| Loading error after sign alignment | 5.55e-11 maximum |
| Cell-score error after sign alignment | 5.82e-10 maximum |
| Center error | 7.55e-15 maximum |
| Native residual | 2.20e-11 maximum |

The independent SciPy 1.18.1 ARPACK calculation uses a sparse centered linear
operator rather than a dense cells by genes matrix. The comparison checks every
component and cell/feature identity, not selected examples. Native scores use
U times singular values; the ATAC LSI output exposes U and standardized U, so
joint analysis must explicitly choose scaling rather than concatenate silently.

The first reference report mistakenly summarized the normalized shared AnnData
buffer as raw counts. The initial raw-count assertion was correct. The corrected
script captures raw total before normalization; its rerun reports 11,786,194.
The earlier report and log remain in the evidence and are superseded only for
this count-summary field.

## Evidence and next step

The archive and manifest retain the native harness/binary, source and library
fingerprints, protocol, native results, independent result, scripts and logs.
The original HDF5 remains the input documented in [MULTIOME.md](../MULTIOME.md);
its SHA-256 is bound in the protocol. Replay scripts retain development paths.

This completes the paired RNA numerical prerequisite for joint neighborhood
analysis. It does not complete multimodal integration. The next comparison
must include RNA-only, ATAC-only and joint neighborhoods, explicit modality
scaling/component selection, independent numerical checking, and measured
biological preservation. The single-donor dataset cannot validate donor
integration, unseen perturbations or clinical predictions.
