# Reproduction boundary

The authoritative driver is `Tools/Omics/Reduction/MetalDistanceBlocks/EndToEnd/run.py`.
It was invoked with the pinned PBMC3K source, PCA plan, current `numivivo-omics`
binary and `scanpy.py`, using `--cohort PBMC3K --cells-expected 2700 --repetitions 3`.
The binary graph and clustering stages used `singlecell-pca-neighbors`,
`singlecell-pca-neighbors-verify`, `singlecell-graph-cluster` and
`singlecell-graph-cluster-verify` with the retained JSON plans. The execution
was on the physical M4 Pro host with the pinned HDF5 library recorded by the
run environment. Large PCA/graph stores remain on that host; this directory
contains the compact reports, plans, receipts and failure control needed to
review the result.
