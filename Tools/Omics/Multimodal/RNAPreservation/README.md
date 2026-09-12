# RNA preservation diagnostic: not promoted

The native weighted graph does **not** pass the prespecified development
no-worsening check against RNA-only neighborhoods. Results differ between
directed neighbor averaging and fuzzy-graph averaging, so neither the positive
secondary result nor numerical correctness establishes biological preservation.

All 2,711 measured nuclei and all 34,601 genes outside the 2,000 selected RNA
PCA genes are evaluated. The source file, complete native product result and
three native baseline graphs are hash-bound before scoring. No gene, cell or
support group was selected after seeing these results.

| Weighted versus RNA-only | Aggregate squared error | Macro variance-normalized error | Gate |
| --- | ---: | ---: | --- |
| Mean of 14 non-self selected neighbors (primary) | +0.8527% | +0.3312% | FAIL |
| Row-normalized fuzzy graph (secondary) | -0.1274% | -0.4175% | PASS |

Positive changes indicate worse reconstruction. The frozen rule required no
increase in aggregate and macro error for both graph uses, plus no increase in
aggregate error within every nonempty detection-frequency group. The primary
comparison worsens in all three groups; the secondary improves in all three.
The combined rule therefore fails. It is a development diagnostic, not a
clinically calibrated noninferiority margin or a statistical significance test.

## Complete comparison

Errors below sum over every evaluated cell/gene coordinate, including measured
zeros. Macro error averages per-gene SSE divided by centered observed sum of
squares for variable genes detected in at least 20 cells.

| Graph | Neighbor-mean SSE | Fuzzy-mean SSE |
| --- | ---: | ---: |
| RNA-only | 6,921,620.49 | 7,062,271.82 |
| ATAC-only | 7,033,918.55 | 7,136,153.22 |
| Equal-weight joint | 6,938,781.24 | 7,055,030.55 |
| Cell-specific weighted | 6,980,640.15 | 7,053,272.96 |
| Leave-one-cell-out global gene mean | 7,264,815.69 | 7,264,815.69 |

The weighted estimates improve aggregate error over the leave-one-out global mean
by 3.91% and 2.91%, respectively. However, their macro variance-normalized errors
are worse than that simple baseline by 2.38% and 3.71%. Aggregate gains must not
be presented as reliable prediction across genes.

There are 10,796 zero-observed genes, 10,797 genes detected in 1–19 cells,
4,114 detected in 20–99 cells and 8,894 detected in at least 100 cells.
The macro comparison includes 13,008 variable genes. Undefined normalized errors
for zero-variance genes are not silently converted into successful predictions.

## Calculation and checks

Counts are read independently from the original 10x HDF5, retaining every RNA
gene in each cell's normalization denominator. Targets are
log1p(10,000 times UMI count divided by complete RNA cell total).
Prediction excludes each query cell. Uniform-neighbor and fuzzy-row operators
are applied in 128-gene sparse blocks; no complete dense cells-by-genes matrix
is constructed.

Every gene's squared error is checked using two algebraic forms, with combined
relative/absolute floating-point tolerances. Maximum absolute disagreement is
1.35e-9. Scalar averaging across all cells for 17 prespecified genes independently
checks the neighbor-mean calculation; maximum disagreement is 4.33e-12.
The full per-gene arrays, all baseline summaries, null comparison, source hashes,
protocol, code, versions and logs are retained in the verified archive/manifest.

The final protocol is `protocol.json`; the earlier transfer-stage snapshots are
retained separately. All final source hashes were verified after transfer and
before scoring. Inputs come from [PairedWorkflow](../PairedWorkflow/README.md)
and [PairedNeighbors](../PairedNeighbors/README.md); restore those archived inputs
or adjust the recorded development paths when replaying.

## What this does and does not establish

These are genes outside the selected PCA set, **not an independent held-out-gene
trial**: HVG selection and normalization used all genes. The data contain one
donor and the same cells used to construct the graph. This measures conditional
RNA reconstruction, not unseen-donor/context perturbation prediction, cell-type
accuracy, ATAC preservation, regulatory causality or clinical outcomes.

The successful source replay and numerical graph qualification remain valid.
Broad biological-preservation promotion remains withheld. Further method changes
must retain this result and be evaluated against separately specified biological
evidence, rather than selecting only the favorable secondary metric.
