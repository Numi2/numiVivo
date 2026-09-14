# Native gzip-wrapped H5AD import qualification

Date: 2026-09-14
Source commit: `edf8407ead6b1af1c0c8037af1c86d8397227943`
Host: Apple M4 Pro, macOS 26.6 (25G72), Swift 6.3.3
HDF5: 2.2.0 (`/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib`)

The fresh native build passed `swift build`. The focused native suite
`SingleCellInterchangeTests` passed all 6 tests. A native `singlecell-h5ad-write`
created a 2-cell × 2-feature CSR AnnData file; Python supplied only the
external gzip wrapper. The fresh native `singlecell-h5ad-import` then decoded
the wrapper, read the HDF5 payload, reconstructed CSR counts `[4, 5, 6]`, and
retained the exact compressed source bytes in `original.h5ad`.

The first smoke invocation recorded the expected missing-library preflight in
the isolated shell. The qualification rerun used the existing host HDF5 library
path above and completed successfully. This receipt qualifies bounded
external-gzip handling for the direct count-projection import path. It does not
qualify preservation of every AnnData slot through the count projection or any
biological outcome prediction.
