# Measured untreated-cell null stress test

The [frozen protocol](PROTOCOL.md) creates ten balanced sham splits of each
qualified untreated population: 1,316 Kang B cells from eight paired donors and
8,581 Hagai fibroblasts from three paired donors. No measured count is changed.
Each split retains all source genes and original donor, batch and library
identities. Its two arms are disjoint subsamples of the same untreated libraries,
not new biological replicates.

This evaluates conditional false-positive behavior under label randomization.
It does not recreate independently sampled donor experiments or assess power
against genuine perturbations. The splits overlap in cells and genes; summary
counts are descriptive, without independence-based confidence intervals or a
general FDR-control claim.

## Completed results, 2026-09-10

All twenty declared splits completed: forty native fits and forty native replays,
plus 120 R fits. There were no failed commands, final native fit nonconvergence,
or R warnings/messages. All forty native count, source-membership, design,
normalization and BH checks passed. The largest independently recomputed Wald
probability error was 4.45e-16; the largest size-factor error was 6.67e-15.
The native annotations preserved all 34 original Kang and ten original Hagai
HDF5 datasets and added the ten exact split columns. Archive restoration passed
native replay; three corrupt-source/payload/schema cases failed before publishing
an output bundle.

The following are method-specific tested families using fixed native factors.
Calls are summed across ten overlapping splits, so they are not unique genes.

| Study | Method | Tested genes per split | Total BH<0.05 calls | Splits with any call |
|---|---|---:|---:|---:|
| Kang | Native default | 3,202–3,204 | 0 | 0/10 |
| Kang | Native active donor | 7,378–7,397 | 0 | 0/10 |
| Kang | edgeR robust QL | 7,424–7,426 | 9 | 3/10 |
| Kang | limma-voom | 7,424–7,426 | 16 | 6/10 |
| Kang | DESeq2 | 7,424–7,426 | 0 | 0/10 |
| Hagai | Native default | 11,917–11,921 | 20 | 8/10 |
| Hagai | Native active donor | 11,917–11,921 | 20 | 8/10 |
| Hagai | edgeR robust QL | 12,078–12,081 | 112 | 10/10 |
| Hagai | limma-voom | 12,078–12,081 | 102 | 10/10 |
| Hagai | DESeq2 | 12,078–12,081 | 1 | 1/10 |

Package normalization gives Kang totals of 9/27/0 and Hagai totals of 111/93/1
for edgeR/voom/DESeq2, respectively. The complete per-seed results retain both
normalization modes, ordinary DESeq2 filtering, withheld genes, raw probability
fractions and separate joint-family comparisons. The native default support
gate's substantial Kang coverage loss remains visible.

**Native Hagai sham calls are observed counterevidence to calibration.** Both
native policies report twenty calls covering fourteen distinct Ensembl genes.
No algorithm or filter was changed in response. Zero Kang calls and fewer calls
than some references do not establish superiority: power, independent donor
sampling and general FDR control were not evaluated here.

`audit_null_hits.py` is explicitly a post-score diagnostic. It independently
reconstructs each affected donor/arm count and retains cell-count concentration,
dispersion, fitted means and Cook's distance. Seven of twenty calls have at least
one nonzero arm with more than half its gene counts from one cell; eight calls
have every arm's largest cell below ten percent. None was flagged as a dispersion
outlier, and maximum Cook's distances span 0.135–15.426. Thus a rare-cell
explanation alone cannot cover every call; these diagnostics do not identify a
causal repair. Robust small-sample inference and independent calibration remain
open, and this benchmark must stay unchanged when evaluating future changes.

The frozen protocol SHA256 is
`4203a2d11df8e58afb973b814055e7d7b8aed128a1aefc5995c785cfc93f1f57`.
The [evidence manifest](evidence/2026-09-10/manifest.json) binds compact inference,
full R tables, exact assignments/count references, logs and checks. Complete
native model reports and original/annotated H5ADs remain at individually hashed
Mac mini paths recorded in that manifest.

## Execution and verification

`prepare.py` freezes SHA256-ranked assignments without examining expression,
adds native annotation plans, and computes only small donor/arm-by-gene reference
matrices. `check_preparation.py` independently reconstructs every assignment and
all twenty sparse count aggregations. `check_annotation.py` verifies every
original HDF5 value, dtype, shape, storage setting and attribute, plus all ten
added columns, against native annotated data. Original condition labels remain.

`run_native.py` uses the existing full release product at
`f40b37c870e6e1702e3feff8f7b400188cce82ca`, executable SHA256
`4eee4da0f3cdfef0ea47d341a6ce75d27042d2ff8e017b57be3f5e9602bb1884`.
Each split runs separately under original full-support and explicit active-donor
profile policies. Native HDF5 streams source-bound selected cells into aggregates.
Every successful publication is reconstructed by native replay before archival.
No trend, support, influence or filtering choice changes after results.

`run_reference.py` reuses the pinned [Bioconductor comparison](../../Bioconductor/README.md)
with robust edgeR QL, robust limma-voom and DESeq2 Wald, each under fixed native
factors and package normalization. DESeq2's ordinary result filtering is retained
separately. Packages, warnings and failures are recorded for every run. The
reference input independently reconstructs paired designs and median-ratio
factors; `summarize.py` verifies them against native results.

`summarize.py` verifies every native aggregate and source-row membership against
the independent reference, all feature identities, designs and normalization,
and independently recomputes BH. It reports method-specific full tested families
and a separate common intersection, keeping withheld genes visible. Explicit
`--allow-incomplete` mode labels missing results as pending; final evaluation
requires every declared split and method. Missing/failed results remain recorded.

## Artifact storage

Successful native bundles become compressed, source-deduplicated archives.
The common annotated H5AD remains retained under each study. Exact source hashes
and open handles are checked before removing duplicate snapshots; gzip round-trips
are checked before removing plain JSON. Each archive records source/payload
hashes and successful native publication/replay. Materialize a fresh normal
bundle before using the CLI verifier:

```sh
python restore_bundle.py --archive /runs/kang/1/default --out /new/native-bundle
numivivo singlecell-h5ad-pseudobulk-verify /new/native-bundle
```

The restore tool verifies the shared source and all payload hashes; existing
outputs are rejected. Detailed native model reports and large source H5ADs may
remain externally retained with hashes. Compact per-gene inference, reference
outputs, assignment plans and summary checks provide reviewable evidence.

## Reproduce

Use the qualified Python/R environments and native HDF5 library. Prepare both
studies before fitting either model family:

```sh
python prepare.py --study kang --source-bundle /qualified/kang --out /runs/kang
python prepare.py --study hagai --source-bundle /qualified/hagai --out /runs/hagai
python check_preparation.py --root /runs --source-root /source-parent --out /runs/preparation-check.json
python run_reference.py --root /runs --r-library /qualified/R-library
mkdir -p /native-runs/kang /native-runs/hagai
for study in kang hagai; do
  cp /qualified/$study/original.h5ad /native-runs/$study/original.h5ad
  cp /runs/$study/native-input.json /runs/$study/annotation-plan.json /native-runs/$study/
done
python run_native.py --root /native-runs --binary /qualified/numivivo --hdf5 /qualified/libhdf5.dylib
rsync -a --exclude='*.h5ad' /native-runs/ /runs/
python summarize.py --root /runs --out /new/final-results
```

`check_preparation.py` uses the existing `numivivo-kang-r-native-current-20260909`
and Hagai counterpart names under `--source-root`; source hashes are checked.
The native execution root must be fresh, with only the three inputs per study;
the runner creates seed directories itself. The reference root already contains
prepared seed directories and is therefore a separate path. Retain shared native
H5ADs in the native root when restoring its archives.
If execution is split across hosts, copy the prepared metadata/source files to
the native host and merge its compressed result archives back before summary.
For this retained run, `summarize.py --remote-models` reads missing detailed model
reports from `macmini:/Users/n/numivivo-null-benchmark-20260910` and verifies both
compressed and logical hashes before using them. This avoids redundant local
model mirrors. `archive.py` packages compact evidence and verifies every stored
and decoded hash; it records external sources instead of implying they are in Git.
Python is authoring/reference tooling, not part of native fitting.
