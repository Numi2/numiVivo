# Complete native HIRISA graph clustering

**PASS for the declared full-cohort operational and numerical experiment.**
The native publication, native replay and independent partition check complete
on all **1,612,594 original cells**. The result has **30 connected communities**.
These are graph partitions, not calibrated cell types or a complete biological
annotation benchmark.

The executable SHA256 is
`67343146e05b783b1bb0806cb765a536890109a6b50816d34617144e0af7c114`,
from the frozen sequential-ridge/runtime build. This is the distinct buffered
reader qualification. The earlier mapped-reader publication has now finished
and produces the identical complete result. Its own reconstruction failed for
lack of disk space; that historical failure remains separate from this passing
buffered-reader qualification.

## Native and independent evidence

| Stage | Elapsed seconds | Result |
| --- | ---: | --- |
| Fresh complete PCA | 615.879 | All six original fitted payloads exactly reproduced |
| Complete neighbor graph | 1160.647 | All seven original graph payloads exactly reproduced |
| Native clustering publication, including parent reconstruction | 1411.683 | Completed |
| Native clustering replay, including parent reconstruction | 1498.159 | Completed |
| Independent every-cell partition check | 9.018 | Passed |

Times were observed on the shared physical Mac mini with concurrent workloads;
they are not controlled method speed comparisons. The experiment preserves its
original [clustering protocol and three igraph references](CLUSTERING_REFERENCE.md).

Every original cell identity and label, every graph edge in the objective and
every community's connectivity were checked. Native modularity is
**0.8857395360221528**, independently recomputed as **0.8857395360221637**.
There are zero disconnected communities. The maximum objective gap to the three
original full igraph partitions is **0.00181198994031051**, below the frozen 0.02
margin. The native run reports 3,089,520,927 edge visits, below its frozen budget.
The objective comparison does not require partitions from different stochastic
methods to have identical labels.

The new graph receipt is
`58d5288c72762c459dae1d645f06bb10b0f681a33fe9606661c09325267f1314`.
Reuse of earlier independent graph/igraph evidence is justified by exact graph
and fitted-PCA payload equality and matching original count-source identity.
The archive retains both old and new receipt identities; it does not relabel
the original references as having run under the newer native executable.

## Evidence and remaining scope

The [full result archive](evidence/2026-09-10-full-clustering/manifest.json)
contains the complete native result with every cell label, receipts, runtime
and phase bindings, publication/replay logs, independent checks and collection
proof. Large byte-identical parent graph/PCA payloads restore through the
manifest's exact paths into the earlier graph/storage archives. The three
original reference partitions remain in the clustering-reference archive.
All 51 stored and decoded members can be checked with:

```sh
python Tools/Omics/Benchmarks/HIRISA/verify_archive.py \
  Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-10-full-clustering
```

This closes this full native clustering/replay gate. Rare-cell preservation,
cluster annotation, prospective reference transfer and phenotype prediction
require separate evidence. In particular, the [program-preservation failures](INTEGRATION_PROGRAMS.md)
remain failures; correct graph arithmetic cannot qualify biological integration.


## Original mapped-reader baseline, completed 2026-09-11

The original executable `ee332d69aa3b3514444f00c4474f03332c296e236ed61d26cf2e2fb933e3d699`
finished publication in **32,819.665 seconds** with peak RSS **4,445,929,472 bytes**.
Its 118,692,184-byte result is exactly identical to the qualified buffered-reader
result, SHA-256 `94c7b1113226676bf865bcebe38939b2447f98b8d56dc054e77d943b7c0b32f3`.
The original plan is also byte-identical. The waiting original coordinator ran
its own complete independent check: all 1,612,594 cell identities/labels, graph
edges in the objective and community connectivity pass, with the same modularity
and three reference comparisons reported above.

The subsequent original reconstruction **failed after 330.650 seconds with
ENOSPC**. The original coordinator and worker have exited; no restart was made.
The historical reconstruction therefore remains failed. The already completed
buffered-reader publication/replay and independently checked identical result
remain valid separate evidence, rather than being used to relabel that failure.

Original publication took 23.249 times the buffered run's observed 1,411.683
seconds. These were different executables on a shared host with changing competing
workloads, not a controlled repeated benchmark. The original log records
31,538.98 seconds of system time and 1,268.74 seconds of user time; the observation
supports investigating storage/mapping overhead but is not isolated causal proof.
Future reconstruction runs need sufficient temporary-disk headroom before launch.

The [terminal baseline archive](evidence/2026-09-11-mapped-baseline/manifest.json)
preserves original publication/failure logs, receipts, process outcomes, independent
checks and source. Its manifest reuses the exact complete result already stored
in the full-clustering archive, avoiding another copy of every label.
