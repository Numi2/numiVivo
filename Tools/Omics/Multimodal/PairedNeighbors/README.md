# Paired RNA, ATAC and joint neighborhood baselines

Native exact neighborhoods and fuzzy connectivity now pass independent
calculation on all 2,711 paired nuclei. This is a measured-data baseline for
learned multimodal integration, not a completed integration method.

The frozen protocol uses all 30 qualified RNA PCA and all 30 ATAC LSI components.
RNA scores are standardized by sample standard deviation; ATAC uses the qualified
standardized U output. Each modality row is then L2-normalized. The joint vector
concatenates the two modality rows divided by sqrt(2), giving equal squared-distance
weight to each modality. No component is silently removed, including the
depth-correlated first ATAC component. Zero or nonfinite modality norms reject.

The actual native graph owner executes through a linked scoped harness with
test visibility. There are 15 neighbors including self, local connectivity 1,
and 3,673,405 unordered distance pairs per graph. Native search retains bounded
heaps, not a dense pair matrix or cells by genes matrix. Independent SciPy cdist
runs in 128-row blocks; UMAP-learn 0.5.12 constructs reference fuzzy connectivity.

| Graph | Non-self directed CSR entries | Components | Isolated cells | Maximum reference weight error |
| --- | ---: | ---: | ---: | ---: |
| RNA | 49,474 | 1 | 0 | 2.71e-6 |
| ATAC | 55,046 | 1 | 0 | 6.63e-6 |
| Equal-weight joint | 52,052 | 1 | 0 | 2.77e-6 |

All neighbor indices match exactly. Maximum distance error is 1.34e-15 and all
fuzzy weights pass the frozen 1e-5 tolerance. Exact paired cell identities and
input hashes are checked. Graphs are symmetric; the entry counts above count
each undirected connection twice.

## Structural diagnostics

| Pair | Mean fraction of 14 shared non-self neighbors | Cells with no shared neighbors |
| --- | ---: | ---: |
| RNA / ATAC | 9.23% | 1,053 |
| RNA / joint | 43.82% | 38 |
| ATAC / joint | 32.13% | 79 |

Equal modality distance weights do not imply equal neighbor retention. These
numbers describe graph structure; they do not establish cell-type accuracy,
biological preservation, regulatory links or useful integration. No inferred
labels are used as ground truth. The single donor cannot qualify batch integration.

## Replay and next work

The archive and per-member manifest retain the frozen protocol, exact native
harness/binary, source-build fingerprint list, all three complete graphs and
scaled representations, independent checker and diagnostics. Inputs are the
published [paired RNA PCA](../PairedRNA/README.md) and [ATAC LSI](../LSI/README.md)
results, bound by SHA-256 in the protocol. Scripts retain development paths;
restore the corresponding inputs before replay. The first harness compile found
a cell versus identity type mismatch; explicit identity conversion fixed it before
any run. The successful compiled harness is retained.

Next, implement and independently compare learned modality weighting against
these three fixed baselines, then test biological preservation using separately
specified experimental evidence. This baseline must not be promoted as learned
weighted-nearest-neighbor integration, million-cell scaling or prediction.
