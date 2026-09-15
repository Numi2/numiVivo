# Current PBMC3K PCA, neighbors and clustering receipt

This receipt records a fresh physical M4 Pro run of the current published
NumiVivo binary (`476ad07914bd6edc1a92b075ca2bbc9057afa035`) on the pinned public
PBMC3K source. The source is the Scanpy PBMC3K file identified by upstream
SHA-256 `89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1`.
The source URL is `https://falexwolf.de/data/pbmc3k_raw.h5ad`.
The native input is an AnnData re-encoding with explicit
`obs["native_sample"] = "pbmc3k-public-library"` and `var["gene_ids"]`; its
source SHA-256 is `d5bed3f71eb771484450369b17e11d4f7df9c38f965d3a8e2b72cda0de28004d`.
It contains 2,700 cells, 32,738 features and 2,286,884 nonzero count values.

The predeclared route ran three fresh repetitions in fixed orders through raw
sparse H5AD → native 2,000-feature Seurat-HVG/20-component PCA → 20-neighbor
fuzzy graph, with independent Scanpy raw-count normalization, HVG, PCA and
neighbors. CPU and Metal each passed PCA/graph verification on every repetition.
The final repetition has 51,300 non-self neighbor entries and 77,916 graph
edges per backend, with zero missing or extra neighbors/edges and exact edge
coordinates. Maximum common distance differences are `1.5122036955972362e-10`
(CPU) and `2.7322685482999987e-06` (Metal); maximum fuzzy-weight differences are
`3.0030269828618117e-06` and `4.541294049920097e-06`. Selected feature IDs,
cell barcodes and sample IDs match in order.

Native seeded binary-graph Louvain replay was then run for seeds 7, 19 and 42.
CPU and Metal produce the same label vector and 100% pair-relation agreement at
each seed, with ten communities and zero disconnected communities. Modularity
is 0.7501332186/0.7501332153 (seed 7), 0.7527419247/0.7527419213 (seed 19) and
0.7477234554/0.7477234536 (seed 42), CPU/Metal respectively. These are
numerical and descriptive graph checks; cluster labels are not cell-type
annotations.

Median native pipeline time is 2.0486928332 s CPU and 2.0197088341 s Metal;
median Scanpy time is 6.1914758747 s. Maximum recorded RSS is 178,061,312,
198,819,840 and 627,015,680 bytes respectively. The timing boundary includes
startup, reads, arithmetic and writes; verification is outside it. Scanpy's
full output digest varies because its timing report is run-specific, while its
deterministic artifacts match across repetitions.

The first attempt is retained in `pre-fix-namespace-comparison.json` and
`scanpy-pre-fix.py`. It failed because Scanpy reported display-name indices while
the native mapping used `gene_ids`, giving zero selected-feature overlap. The
explicit namespace selection in the current `scanpy.py` fixes that contract; the
failed evidence remains part of the audit.

This receipt establishes current numerical interoperability and descriptive
clustering for one unreported PBMC3K library. It does not establish donor
replication, biological preservation, perturbation prediction, phenotype or
clinical outcome prediction.
