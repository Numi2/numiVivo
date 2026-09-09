# File-backed multilevel Louvain

The binary PCA graph store now feeds the same deterministic Louvain solver as
the resident count-analysis route. Original edges and every aggregated level
use CSR files with 16 MiB reader windows. Cell identities, offsets, labels,
degrees, community totals, permutations, traversal queues and output JSON remain
resident. A single aggregate row's edge map also remains resident; a high-degree
row can therefore still use memory proportional to the number of communities.

```sh
numivivo singlecell-graph-cluster /path/to/binary-graph \
  --plan clustering-plan.json --output /new/clustering-bundle
numivivo singlecell-graph-cluster-verify /new/clustering-bundle
```

```json
{
  "schemaVersion": 1,
  "clustering": {
    "resolution": 1,
    "seed": 7,
    "maximumSweeps": 100,
    "maximumLevels": 32,
    "levelTolerance": 0.0000001
  },
  "maximumEdgeVisits": 1000000000
}
```

The input must use `storage: "binary"`. Fitted and frozen-query PCA graphs,
from either exact or HNSW search, are supported. Publication retains and
reconstructs the complete parent graph and PCA lineage before clustering.
Verification reconstructs the lineage and partition, then checks all result
hashes. Receipts bind the executing binary; a new executable must refit inputs
from their original sources rather than reuse an older executable's receipts.

Aggregation first counts records per new community, then scatters source edges
into a temporary file in their original row/column order. Each new row is reduced
separately and written in sorted column order. This preserves floating-point
addition order and both directions of internal edge mass on the diagonal.
No global edge array or dictionary of all aggregate edges is constructed in
file mode. Private level files are removed on success and failure.

`maximumEdgeVisits` defaults to one billion and permits at most 20 billion.
It includes repeated graph reads and temporary aggregation-record visits;
it is separate from the parent PCA and neighbor-search work budgets. Cancellation,
work exhaustion, failure to converge and objective decrease prevent publication.
Input axes permit at most one million cells and 254 million directed CSR records;
these admission limits are not million-cell qualification.

`check_file_clustering.py` checks fitted/query inputs in exact/HNSW modes,
byte-identical replay, a frozen previous-algorithm result, and rejection of
rehashed labels, execution counters and parent edges, insufficient work budget,
JSON input and overwrite. The frozen oracle is built only by
`build_legacy_clustering.sh`; it is never linked into the product. Its metadata
records the exact pre-refactor source revision and transformation.

The full-cohort protocol requires every input graph coordinate and FP64 bit to
equal the previously independently qualified graph, and the entire result JSON
to equal the frozen solver. NetworkX independently evaluates modularity and
connectivity and runs Louvain with seeds 7, 19 and 41. The predeclared maximum
objective deficit is 0.02. Different heuristic partitions can pass this numerical
comparison. Metadata associations remain descriptive, with absent or single-level
annotations reported explicitly.

[File-backed embedding](FILE_EMBEDDING.md) now shares the binary graph input.
File-backed integration, million-cell execution, partition
stability across resolutions, stronger biological evaluation, Leiden refinement
and Metal/scverse end-to-end comparisons remain open.

## Full-cohort qualification

The [evidence archive](evidence/2026-09-09-file-clustering/README.md) retains
all full-cohort graph checks, result bytes, execution counters, three reference
partitions per cohort, source/binary fingerprints and command logs. All three
complete result JSON payloads equal the frozen pre-refactor algorithm exactly.
Every graph coordinate, FP64 value, PCA score byte and cell identity also equals
the previous qualified input. Norman retains the earlier HNSW recall scope of
2,048 exact query comparisons; that sampled search-quality result is separate
from the complete graph/partition checks here.

| Cohort | Cells | Clusters | Native modularity | Best reference deficit | Native/reference ARI range |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baron | 8569 | 19 | 0.872844496 | 0.000339582 | 0.867–0.949 |
| Hagai | 13863 | 11 | 0.800139931 | 0.002965046 | 0.711–0.764 |
| Norman | 111445 | 17 | 0.748492374 | -0.004046131 | 0.431–0.486 |

Independent objectives agree within 3.4e-14. Disconnected-community counts
are 0, 0, 0 for Baron, Hagai and Norman. All reference deficits pass the
predeclared 0.02 threshold. This does not establish equivalent partitions or
authoritative cell types. Baron cell-type ARI is descriptive (0.496); Hagai lacks
cell-type and donor annotations in this prepared input, and Norman has a single
K562 group and no declared independent donors.

| Cohort | Publish seconds | Publish peak resident bytes | Verify seconds | Verify peak resident bytes |
| --- | ---: | ---: | ---: | ---: |
| Baron | 7.38 | 205389824 | 7.39 | 207683584 |
| Hagai | 16.19 | 224493568 | 16.07 | 221118464 |
| Norman | 264.95 | 482918400 | 261.68 | 475971584 |

These are single-run native CPU observations on M4 Pro/24 GiB/macOS 26.6.
Each publish/verify reconstructs the entire parent PCA and graph; the isolated
frozen oracle parses an already-built JSON graph and its timings are not a
workflow speed comparison. Independent NetworkX checks ran on the laptop.
Completed artifacts were copied while subsequent native cases ran; no native
build or competing benchmark ran concurrently. No same-host scverse or GPU
speedup claim is made.

The release and scoped builds passed, as did 80 Swift tests in 24 suites,
covering SingleCell, MultiAssay and CountStore, 12 full-cohort production
commands, three frozen-oracle executions and 66 CLI
fixture commands (including 25 expected rejections). The archive retains the
initial test compilation failure and the successful correction. Numerical
qualification remains scoped to the recorded binary and production-source hashes.
