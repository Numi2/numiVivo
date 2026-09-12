# Native source-to-graph paired multiome workflow

The product CLI now runs the complete paired workflow from an original 10x HDF5
count file: native import, RNA normalization/HVG/PCA, ATAC TF-IDF/LSI,
cell-specific modality weights, directed weighted neighbors and fuzzy connectivity.
It retains the original source, explicit plan, complete results and receipt.
A separate verifier reconstructs the complete result from the retained source.

```sh
bash Tools/Omics/H5AD/build.sh /tmp/numivivo-omics --with-cli
/tmp/numivivo-omics/numivivo-omics multiassay-10x-paired original.h5 \
  --plan Tools/Omics/Multimodal/PairedWorkflow/paired-plan.json --output paired-bundle
/tmp/numivivo-omics/numivivo-omics multiassay-paired-verify paired-bundle
```

A native HDF5 library must be discoverable, or set `NUMIVIVO_HDF5_LIBRARY` to
its absolute path. Qualification used the installed HDF5 2.0.0 dylib bundled with
h5py 3.16.0; the exact path and library hash are recorded in `runtime.json`.
Python did not perform the product's import, transforms or graph calculation.

## Plan and behavior

The supplied plan maps the complete [10x paired dataset](../MULTIOME.md), with
RNA UMI counts and ATAC cut-site counts kept distinct. It requests 2,000 RNA HVGs,
30 RNA PCs, 30 ATAC LSI components and 15 neighbor slots including self.
Normalization remains fixed at 10,000 over all RNA genes; ATAC uses the qualified
count-based method-1 TF-IDF. Each modality is standardized and row-L2-normalized
before the qualified weighting/graph calculation. No component is silently
excluded, including the depth-correlated first ATAC component.

The plan rejects unknown keys, invalid schemas, incompatible assay mappings and
invalid reduction/neighbor options. Reduction options are validated before their
component counts are combined, preventing overflow from extreme plan values.
RNA and ATAC must have exactly matching observation-index arrays. If default
positive-count/feature RNA filtering would remove any paired cell, the workflow
rejects rather than silently changing the pairing. No mitochondrial filtering
or inferred cell annotation is performed.

The bundle contains:

- `original.h5`: immutable source snapshot.
- `plan.json`: explicit mapping and algorithm settings.
- `result.json`: RNA metadata and reduction, ATAC metadata, full LSI result and
  count statistics, modality diagnostics, and weighted graph.
- `receipt.json`: hashes binding source, plan, result and producing executable.

Publication uses a private staging directory and occurs only after success.
Existing output destinations are preserved. Verification reports its own
executable fingerprint separately from producer provenance and compares the
reconstructed result bytes, not only user-supplied checksums.

## Complete measured qualification

All 2,711 nuclei, 36,601 RNA genes and 98,319 peaks are retained.
The native product reproduces every previously qualified RNA PCA field, ATAC
loading/U/standardized embedding, ATAC axis/statistic field and weighted-graph
field exactly. Inputs include 5,218,473 RNA nonzeros / 11,786,194 UMIs and
19,292,713 ATAC nonzeros / 48,245,242 cut sites.

The complete source replay passes. Rejection/preservation checks cover unknown
keys, schema mismatch, incompatible assay selection, overflowing component
requests, existing destinations, and a valid-JSON result whose method text and
receipt checksum were both modified. Source reconstruction rejects the latter.

The initial local compile passed; a later final compile exhausted local temporary
storage. The final build and complete tests ran in an isolated Mac mini directory.
Its first execution lacked a configured HDF5 path and rejected before publication;
the successful run explicitly selected the installed native library. Both failure
logs are retained, alongside executable, compiler, runtime and input fingerprints.

The archive and per-member manifest retain complete output, receipt, tests,
actual source and build/runtime evidence. Original HDF5 and prior reference
outputs are hash-bound and available through the preceding benchmark evidence;
they are not duplicated in this archive. Replay scripts retain development paths.
Some local duplicates were removed only after published-archive byte verification;
restore archived inputs when replaying older harnesses.

## Qualification limits and next target

This workflow is bounded and resident: the native multiassay reader retains sparse
counts, TF-IDF retains one value per ATAC nonzero, and reductions/results remain
resident. The source reader's 32-million-nonzero limit, exact neighbor-pair budget
and 512-MiB result limit apply. This does not qualify million-cell execution,
peak memory, GPU execution or comparative speed.

Full Seurat approximate-search/SNN-graph parity is not claimed; the product uses
the previously qualified exact weighted-neighbor/fuzzy-union method. Cell labels
are not authoritative and this single donor cannot validate donor integration,
unseen perturbations, regulatory links or clinical outcomes.

The [RNA preservation diagnostic](../RNAPreservation/README.md) now fails its
combined no-worsening rule: primary neighbor averaging worsens versus RNA-only,
while fuzzy averaging improves slightly. Numerical/source replay remains valid,
but biological promotion is withheld. Independently specified experimental
evidence is still required for reliable biological outcome prediction.

The subsequent [ATAC diagnostic](../ATACPreservation/README.md) also fails its
combined gate. Neither modality currently supports preservation promotion.
