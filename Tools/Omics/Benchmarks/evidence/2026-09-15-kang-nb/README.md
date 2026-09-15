# Kang 2018 paired-donor NB2 comparison (2026-09-15)

The current native CLI ran the real Kang 2018 B-cell source through sparse
preparation, count import, pseudobulk aggregation, QC, normalization, native
negative-binomial DE and replay. PyDESeq2 received the identical 16 donor-condition
pseudobulks and the same predeclared count filter. No cells or genes were
reconstructed from normalized expression, and cells were not treated as
independent biological replicates.

The run passed for 2,651 cells, 15,706 source genes,
1,479,543 nonzeros, eight donors and 16 pseudobulks. Counts,
pseudobulks and Scanpy QC were exact; the maximum normalization error was
8.88e-16. Native NB2 tested
5,400/8,894 eligible genes;
3,494 were retained as explicit
rank-deficient support rejections. On the common tested set, effect Spearman
correlation was 0.999067, sign agreement was
99.00%, and the top-50 BH overlap was
43. All five predeclared IFNB response genes were
positive in both models.

`report.json` retains the per-feature diagnostics and reference comparison;
`execution.json` is the compact receipt; `sources.sha256` binds the scoped
native source set. The downloaded H5AD, native store and large result tables
are intentionally not committed. This is one predefined real paired-donor
contrast and a descriptive model comparison. It does not establish FDR or
interval calibration, cross-study transfer, causal mechanism, or general
biological-outcome prediction.

Reproduce from the repository root with the pinned benchmark environment:

```sh
python Tools/Omics/Benchmarks/run_kang.py \
  --binary /tmp/numivivo-kang-build/numivivo-omics \
  --source /tmp/numivivo-kang.h5ad \
  --out /tmp/numivivo-kang-run \
  --model negativeBinomial
```

The source is `https://ndownloader.figshare.com/files/34464122` with SHA-256
`e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830`. Versions: Scanpy 1.12.4,
AnnData 0.13.3.post0, PyDESeq2 0.5.4,
NumPy 2.5.3, SciPy 1.18.1,
pandas 3.0.5, h5py 3.16.0;
native HDF5 2.2.0.
