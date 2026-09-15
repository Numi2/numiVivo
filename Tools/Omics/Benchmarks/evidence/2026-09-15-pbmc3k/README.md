# Public PBMC3K native route (2026-09-15)

The current repository revision ran the actual public PBMC3K H5AD through
annotation, native count import, count processing, analysis and replay, then
compared reconstructed counts, QC totals and log-normalization values with
Scanpy. No cells or features were selected or dropped.

The run passed for 2,700 cells, 32,738 features and
2,286,884 nonzeros. Counts and Scanpy QC were exact; the maximum
normalization error was 8.88e-16. Peak RSS was
1,413,054,464 bytes during replay.

`report.json` is the compact result receipt and `sources.sha256` binds the
scoped native source set. The downloaded H5AD, native artifact store and logs
are intentionally not committed. This is a real public-file interchange and
count/QC benchmark, not a donor-replicated biological prediction benchmark.

Reproduce from the repository root with the pinned benchmark environment:

```sh
python Tools/Omics/H5AD/check_public_data.py \
  --binary /tmp/numivivo-h5ad-build-pbmc/numivivo-omics \
  --out /tmp/numivivo-pbmc-run \
  --full-product --count-analysis
```

The script downloads and verifies `https://falexwolf.de/data/pbmc3k_raw.h5ad` with SHA-256
`89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1`. The receipt records AnnData 0.13.3.post0,
h5py 3.16.0, Scanpy 1.12.4 and native HDF5 2.2.0.
