# Distance-tail exact refinement: rare-cell recovery, program gate fails

The hybrid matcher recovers Kang's per-type recall and condition gates, but
**fails within-stratum program preservation**. Spearman correlation is
**0.535757**, versus original PCA's **0.587057**: the 0.051300 loss exceeds the
unchanged 0.05 allowance. Hagai and Ding pass their measured gates. No method is
promoted, and no threshold is relaxed.

| Complete original cohort | Cells | Exact-anchor recall | Refined queries | Fit seconds |
| --- | ---: | ---: | ---: | ---: |
| Hagai | 13,863 | 100.000% | 978 | 3.839 |
| Kang | 24,673 | 99.993% | 2,145 | 15.893 |
| Ding | 44,031 | 99.939% | 3,705 | 33.167 |

All **82,567 cells** remain in fitting and evaluation. A query means a cell plus
its lower- or upper-level search direction. Kang megakaryocyte recall improves
to **0.778153** against baseline 0.750000, and condition balanced accuracy is
0.974326 against baseline 0.973819. These passes do not cancel the program loss.
Four original Kang classifier strata remain unavailable; Ding labels remain
partial. Earlier approximate-matching failures remain retained.

## Fixed, label-free intervention

Start with the width-512 eligible-prefix/suffix HNSW candidate. Within each
observed donor/batch level and eligible search direction, select the ceiling of
5% of queries with greatest approximate kth-neighbor Manhattan distance. Break
selection ties by ascending original row ID. Recompute those queries exhaustively
against every eligible row using grouped-source-index distance ties; leave every
other query unchanged. M=16, construction ef=200, seed=7, k=20, panorama assembly
and the qualified full-anchor Gaussian sigma=15 kernel stay fixed.

No labels, reference anchors or outcome measurements enter this selection rule.
Its design follows inspected earlier failures, so this is development reuse,
not an untouched validation. Large neighbor distance is a heuristic for sparse
local geometry, not an exactness certificate for unrefined rows. The native
refinement has an explicit distance budget and leaves outputs invalid on an error.
The Python research wrapper records selected row IDs and actual refinement work.
This is not a production Swift API or million-cell qualification.

## Execution and independent checks

Actual native fitting and full replay pass on all three complete cohorts.
An independent SciPy cityblock/lexicographic checker verifies all **6,828** selected
queries' exact IDs and distances, with zero distance discrepancy. It reproduces
the selection rule and confirms all unselected outputs are unchanged. Every
selected matching distance and every full-Gaussian update also pass the existing
independent checker. Focused native tests cover ties, level eligibility, exhausted
budgets, insufficient output capacity and nonfinite inputs.

The initial refinement checker repeatedly decompressed NPZ arrays inside its
query loop. That verification attempt was deliberately interrupted and retained;
the replacement loads each array once with unchanged numerical operations and
acceptance rules. Fitting was not restarted. The archive preserves both checker
sources, the interrupted status/reason and the completed verification. Timing
values above exclude witness serialization/hashing; they do not establish an
end-to-end application or scverse speedup.

## Evidence and reproduction

`python3 verify.py` checks all archive hashes, fit-time source and library
identities, output/array bindings, complete cohort accounting, native tests and
the independent refinement receipt. [summary.json](summary.json) retains every
gate. Source, native bridge/header/library, protocol, commands, logs and all JSON
results are archived. [manifest.json](manifest.json) binds full NPZ witnesses
retained remotely, including every selected-query and Gaussian-step array.
Original cohort artifacts, evaluation inputs and qualified Gaussian library
remain separate dependencies; this archive is not a standalone raw-data replay.

Restore those dependencies and run archived `FullGaussianMNN/run.py` with the
recorded library, source manifest, `PROTOCOL_HYBRID.md` and a fresh output path.
Then run `check_refinement.py`, `FullGaussianMNN/check_numerics.py`, and
`LocalMNN/evaluate.py` with the recorded spec. Adjust absolute host paths while
preserving input hashes. Retain the source snapshot under `fit-source`.

The complete biological gate remains false. Improving rare-cell matching while
losing program structure is insufficient; further tuning on these cohorts would
remain development reuse and require fresh independent validation.
