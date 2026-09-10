# Local MNN development result: faster large-cohort fits, failed preservation

The fixed HNSW/64-anchor local correction completed every original Hagai, Kang
and Ding cell and replayed exactly. It retains over 99% of the exact mutual
anchors, but **fails biological-preservation margins on Kang and Ding**. This
candidate is not promoted to the production integration API or defaults.

This is development evidence on already inspected cohorts. It does not establish
independent biological validation, prospective outcome prediction or full HIRISA
integration. Original counts, PCs and reference results remain unchanged.

## Complete-cohort results

| Cohort | Original cells | Exact-anchor precision / recall | Candidate fit / replay seconds | Biological preservation |
| --- | ---: | ---: | ---: | --- |
| Hagai | 13,863 | 99.896% / 99.862% | 5.103 / 5.087 | All applicable original gates pass |
| Kang | 24,673 | 99.994% / 99.948% | 16.325 / 16.318 | Type accuracy, per-type recall and within-stratum program fail; four classifier strata remain unavailable |
| Ding | 44,031 | 99.748% / 99.014% | 26.787 / 26.707 | Type accuracy, per-type recall, global and within-stratum program gates fail; source labels remain partial |

All three improve the original donor/method mixing measure. That improvement
cannot substitute for biological preservation. Kang type balanced accuracy falls
from 0.910778 to 0.876995, beyond the fixed 0.02 allowance. CD8 T-cell recall drops
0.095645 and megakaryocyte recall 0.134291, beyond the 0.05 per-type allowance.
Within-stratum program Spearman falls from 0.587057 to 0.478409.

Ding type balanced accuracy falls from 0.660902 to 0.624435. Recall losses exceed
0.05 for CD14+ monocyte (0.050058), cytotoxic T cell (0.105672), megakaryocyte
(0.057998) and plasmacytoid dendritic cell (0.122629). NK-associated global program
Spearman falls from 0.727308 to 0.674756; within-stratum NK-associated and T-receptor
Spearman fall from 0.398147 to 0.299870 and from 0.201722 to 0.144168. All original
unavailable strata and Ding's 29,411 assigned labels among 44,031 cells are retained.

The [earlier exact owner timings](../MNN_TREE.md) were 4.509, 31.124 and 68.607
seconds respectively. The development driver is slower on Hagai and faster on
Kang/Ding in these observations. These are different driver implementations on
the same physical M4 Pro, not repeated end-to-end CLI speedup measurements.
Candidate timing excludes writing/hashing diagnostic witnesses. The entire
three-cohort fit/replay command took 114.55 seconds with peak RSS 325,287,936 bytes.
The independent original HIRISA clustering baseline remained active.

## Method and retained failed attempt

The [first protocol](PROTOCOL.md) used a single HNSW index with lower/upper-level
result filters. Hagai and Kang fitted and replayed, but Ding exhausted the fixed
500-million-distance-evaluation index budget. Its source, completed outputs,
incomplete Ding state and failure are retained. No biological scores were read
before the [eligible-level follow-up](PROTOCOL_V2.md) was declared.

The follow-up builds an ascending prefix index and a descending suffix index,
querying each level before adding its cells. Each index contains only eligible
reference levels. HNSW remains fixed at M=16, construction ef=200, query ef=128,
seed=7 and 500 million metric evaluations per index lifetime. Matching uses native
FP64 Manhattan distance on normalized original PCs. Each correction step indexes
all source-anchor records, preserving duplicates, then averages biases using the
Gaussian weights of up to 64 returned nearby records per target cell. All cells
participate; this is a local kernel, not a certified approximation to the original
all-anchor Gaussian sum. All other numerical and biological settings stay fixed.

Every candidate alignment order equals its original exact order. Nevertheless,
full-coordinate relative Frobenius differences from exact MNN are 0.174910,
0.251986 and 0.522849. Near-exact anchors alone therefore do not establish retained
coordinates or biological signal. This combined experiment does not isolate the
causal contributions of changed anchor membership and local smoothing.

## Checks, evidence and reproduction

All output-array hashes and fitting reports replay exactly, excluding timing
fields. An independent checker reconstructs every recorded step across 22,081,
237,060 and 391,107 target-row updates respectively. Every selected Manhattan and
squared-Euclidean distance matches exactly; maximum independent delta discrepancy
is 4.45e-16, and final coordinate reconstruction is exact. Global anchor precision
and recall exceed the predeclared 0.95 numerical gate in all cohorts. That gate is
separate from the failed biological gates.

Small exhaustive tests cover both native metrics, filtered/eligible-level
queries, replay and budget/nonfinite rejection. The initial global-index tests pass locally and on the Mac mini; the
eligible-level implementation passes on the Mac mini. Local address/undefined-behavior sanitizers also pass. The reference
metadata, original PCA bytes, protocols, native library and vendored hnswlib
sources are hash-bound before fitting or evaluation.

The [archive](../evidence/2026-09-11-local-mnn) retains both attempts, corrected
scores, mutual anchors, matching neighbors, full biological reports and neighbor/
prediction arrays, numerical checks and source. Bulky per-step source/bias/query/
selection/update witnesses remain externally retained with exact study-relative
paths and hashes. They were checked in full; they are not sampled diagnostics.
The manifest records original input identities and all external restoration needs.

The archive has 68 members and 70,567,851 stored bytes; manifest SHA-256
`76390093170fb4067434651f8c8eb8b860e73ccaf65a406edbcb7673c636b15a`.
After archiving, 53 full step-witness files were relocated to the laptop study
root. The [storage receipt](../evidence/2026-09-11-local-mnn-storage.json) records
all paths and hashes, verification before and after removal, and the 586,657,792
byte observed increase in Mac mini free space. No evidence was discarded; restore
those files from the receipt's destination root when reproducing on another host.
The archive predates this relocation, so use this receipt for current placement.

Build this research driver from a checkout with the recorded vendored sources:

```sh
xcrun clang++ -std=c++23 -O3 -shared -fPIC -I Sources/NumiVivoCore \
  Tools/Omics/Reduction/LocalMNN/LocalNeighbors.cpp -o LIBRARY.dylib
python Tools/Omics/Reduction/LocalMNN/test_neighbors.py --library LIBRARY.dylib --out TESTS.json
OPENBLAS_NUM_THREADS=1 python Tools/Omics/Reduction/LocalMNN/fit.py \
  --source ORIGINAL_COHORT_ROOT --manifest ORIGINAL_ARTIFACTS.json \
  --library LIBRARY.dylib --protocol Tools/Omics/Reduction/LocalMNN/PROTOCOL_V2.md --out NEW_FIT
python Tools/Omics/Reduction/LocalMNN/check_numerics.py \
  --source ORIGINAL_COHORT_ROOT --root NEW_FIT --out NUMERICAL.json
python Tools/Omics/Reduction/LocalMNN/evaluate.py \
  --source ORIGINAL_COHORT_ROOT --root NEW_FIT --spec EVALUATION_SPEC.json --out NEW_EVALUATION
```

Use the frozen source cohort files, original evaluation NPZs and reference checks
listed by the manifest and evaluation spec; adapt absolute paths without changing
identities. Native fitting uses NumPy and the compiled bridge; biological
assessment uses the existing pinned SciPy/scikit-learn evaluators. Do not rerun a
completed fit because an observation expires, overwrite an output directory,
retune these settings or report this candidate as qualified.

Next work should preserve the all-anchor smoothing contract or quantify the
error of an alternative, while retaining the unchanged biological gates. The
current local-kernel failure gives no basis to promote its faster large-cohort
behavior or to skip the remaining complete HIRISA and independent-context tests.
