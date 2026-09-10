# Crowell cortex: eight independent animals

This third treatment study adds an independent-animal design to the earlier
paired Kang and Hagai comparisons. Seven populations complete native analysis,
native replay and 42 direct R reference fits. Exact source/Scanpy count and QC
checks pass. Eight strict conditional-refit p-value checks remain failed, with
their numerical explanation retained below. This is experimental comparison
evidence, not statistical calibration or production qualification.

## Source and frozen scope

The [Crowell et al. paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7705760/)
and [primary sample deposit](https://doi.org/10.6084/m9.figshare.8976473.v1)
identify four Vehicle and four LPS mice, with cortex nuclei and author cell-type
annotations. The muscData `Crowell19_4vs4.Rda`
[object](https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/muscData/Crowell19_4vs4.Rda)
is already filtered by its authors. We retain all **25,224 deposited post-QC
nuclei, 11,076 genes, 38,142,743 nonzeros and 67,349,705 UMIs**.
No additional cell or gene subset is used for import/QC.

`export_counts.R` selects the integer count assay and every row/column metadata
field. `prepare.py` writes a sparse cells-by-genes AnnData matrix and verifies
the complete count round trip, 32 column metadata fields and two row fields.
Character metadata are compared after explicit pandas string representation
normalization; factor levels/order and numeric/logical values remain exact.
This is a count-assay AnnData view, **not a complete SCE conversion**: the
original Rda retains `logcounts` and PCA/TSNE/UMAP. It is also not the unfiltered
Cell Ranger release. Full source Rda and H5AD remain external and hash-bound.

Rda SHA256: `3193a33fa7acb650fdb1e682822705d7c68057aba93d84dac63a364e40fdb1fe`.
H5AD SHA256: `a7ffac5f541409c4b0b55446f6b46f9c955e370b03aa2bc1f785be21c1353b66`.

The deposited workbook independently verifies Vehicle animals LC016/019/022/025
and LPS animals LC017/020/023/026. No pairing or batch is inferred. The design is
intercept plus LPS-minus-Vehicle, with at least ten nuclei per animal/population,
three animals per arm, ten total gene counts and expression in three pseudobulks.
Native median-ratio offsets, Gamma parametric dispersion trend, default support
policy and a fixed normal effect prior of one log2 unit were frozen before fits.
There are 63 observed animal-by-population pseudobulks; LC016 has no CPE nuclei.
CPE retains only two animals per arm after the cell gate and correctly fails
publication with insufficient biological replication. All other populations
retain four animals per arm and six residual degrees of freedom.

## Results

Native execution used the physical M4 Pro Mac mini, CPU FP64, source commit
`3823fb121d41ccd82efb0683bb4a787bde662a90`, and the previously qualified binary
SHA256 `2d9e524cbb1cdd552fc9f52a1c5237cc904f9faec6eb6bc30b7f49d479ee5433`.
The native runtime is unchanged by this benchmark. Seven publications and seven
replays pass; the eighth publication retains the expected replication failure.
Every native report's complete cell QC, identities, all pseudobulk counts and
normalization factors were checked independently; maximum size-factor error is
5.55e-16. No dense cells-by-genes count matrix is constructed.

The table uses fixed native normalization for all methods. Counts are BH<0.05
within each method's own tested family, not comparable error-rate estimates.
Support-rejected native genes remain in the reference family; the complete
comparisons also retain intersection-recomputed BH and unavailable-gene calls.

| Population | Nuclei | Native tested | Support rejected | Native calls | edgeR QL | limma-voom | DESeq2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Astrocytes | 1,740 | 10,650 | 25 | 463 | 509 | 512 | 367 |
| Endothelial | 743 | 10,140 | 22 | 474 | 613 | 548 | 415 |
| Microglia | 416 | 9,435 | 35 | 150 | 297 | 191 | 162 |
| Oligodendrocytes | 1,492 | 10,387 | 7 | 97 | 64 | 64 | 72 |
| OPC | 729 | 10,074 | 16 | 97 | 169 | 142 | 69 |
| CPE cells | 100 | unavailable | — | — | — | — | — |
| Excitatory neurons | 16,622 | 11,070 | 5 | 210 | 27 | 7 | 207 |
| Inhibitory neurons | 3,382 | 10,983 | 10 | 15 | 20 | 6 | 21 |

All 42 R fits complete without captured warnings/messages or DESeq2 coefficient
nonconvergence. Both fixed native and package normalization are retained.
R 4.6.1 uses edgeR 4.10.5, limma 3.68.5 and DESeq2 1.52.0. Scanpy 1.12.4 and
AnnData 0.13.3.post0 supply preparation and independent checks. Full session
information is archived. DESeq2 default filtering/Cook policies are reported
separately from its unfiltered common-design comparison.

Across populations, native effect Spearman correlation with fixed-normalization
DESeq2 is 0.99874–0.999998, with edgeR 0.99418–0.99950 and with limma
0.94346–0.96692. These effect agreements do not resolve the substantial
significance differences, especially excitatory neurons. The predeclared
Stat1/Irf7/Isg15/B2m/H2-D1 panel is positive in all jointly testable native and
fixed-normalization DESeq2 comparisons; Isg15 and B2m are support-rejected in
microglia, and Isg15 in OPC. Positive effects often lack significance. The panel
is descriptive and does not establish biological truth, power or recovery of
the paper's full analysis.

## Retained strict-check failure and diagnostic

`check_models.py` independently fits **72,739 gene/population pairs**, each with
an unpenalized coefficient fit and conditional count-likelihood MAP fit:
**145,478 conditional fits**. Final dispersions are held fixed, so these checks
cannot qualify dispersion estimation or posterior coverage. Coefficients, means,
information-based uncertainty, likelihood and scaled scores pass their declared
bounds. Independent BH arithmetic passes for all populations. All original
SciPy solves report success.

Eight unpenalized refits exceed the predeclared absolute p-value difference
threshold of 1e-7 (six endothelial, one OPC, one inhibitory-neuron gene). Maximum
difference is 1.83177e-7. **The original checker exits 1 and its summary remains
`failed`.** No native optimizer or checker threshold was changed to erase it.

`diagnose_probability.py` separates stored-parameter Wald arithmetic from refit
accuracy. Stored p-values match SciPy to 3.33e-16. A first-order Newton correction
explains the eight refit differences to 2.32e-14, consistent with the existing
1e-7 scaled-score stopping rule. Recomputed BH values differ by at most
2.42e-7 and change no BH<0.05 decisions. This diagnosis is not a replacement
passing result for the strict gate. A diagnostic HYBR progress flag for
ENSMUSG00000036941 is retained: its independent scaled score is 1.45e-16.
The first diagnostic attempt stopped on that flag; the final diagnostic records
both the flag and a separately enforced 1e-9 stationarity requirement.

Preparation also retains two failed attempts: singleton factor levels serialized
as a scalar, and pandas character dtype reconstruction. Both adapter issues were
fixed before native execution. Logs and exact attempted scripts are retained.

## Reproduction and evidence

[Evidence](evidence/2026-09-10) includes frozen plans/protocol, exact pseudobulk
inputs, native feature tables, all R tables, per-gene numerical checks, comparison
summaries, receipts and compressed logs. `manifest.json` verifies stored and
decompressed hashes. Full source matrices and full native reports remain at the
recorded external paths; their hashes and remote report paths are inventoried.
Committed feature tables do not substitute for those complete reports.

With the original Rda and metadata workbook placed in `RUN/source`, use Python
3.12 with the versions above and the R library versions in the session receipts:

```sh
Rscript Tools/Omics/Benchmarks/Crowell/export_counts.R RUN/source/Crowell19_4vs4.Rda RUN/wire
python Tools/Omics/Benchmarks/Crowell/prepare.py --root RUN
python Tools/Omics/Benchmarks/Crowell/run_native.py --root RUN --binary NATIVE_BINARY --hdf5 HDF5_DYLIB
python Tools/Omics/Benchmarks/Crowell/prepare_reference.py --root RUN
python Tools/Omics/Benchmarks/Crowell/run_references.py --root RUN --bioconductor Tools/Omics/Bioconductor --rscript RSCRIPT --r_library R_LIBRARY --python PYTHON
python Tools/Omics/Benchmarks/Crowell/check_models.py --root RUN
# The frozen strict check above exits 1; preserve it and run the diagnostic separately.
python Tools/Omics/Benchmarks/Crowell/diagnose_probability.py --root RUN
python Tools/Omics/Benchmarks/Crowell/summarize.py --root RUN
python Tools/Omics/Benchmarks/Crowell/archive.py --root RUN --out NEW_EVIDENCE
python Tools/Omics/Benchmarks/Crowell/archive.py --out NEW_EVIDENCE --verify
```

Preparation requires a fresh output root beyond its source/wire directories.
The native runner deliberately pins the qualified executable hash. For original
artifact replay, obtain each inventoried external `report.json.gz` alongside
the committed native plan/receipt/archive files and exact `counts.h5ad`, then
use `Tools/Omics/NegativeBinomial/EffectShrinkage/restore_bundle.py` before
`singlecell-h5ad-pseudobulk-verify`.

Runtime and peak-memory logs describe their individual command boundaries.
Native end-to-end publication times must not be compared with R fit-only times
as a speedup. No Metal acceleration was exercised. Cross-study FDR calibration,
the existing Hagai null concern, learned effect priors, posterior coverage and
the broader single-cell objective remain open.
