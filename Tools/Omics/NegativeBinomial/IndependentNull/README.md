# Independent-donor null benchmark

This benchmark uses the complete public Human Immune Health Atlas B/Plasma
release: 160,632 cells, 108 original donors, 32,357 raw features and 233,853,565
stored entries. The 1,199,390,735-byte original H5AD is retained on both hosts.
Every float-stored raw count is an exact nonnegative integer; none is rounded.
The complete data contain 774,987,909 UMIs and 1,772,451 nonzero donor aggregates.

[The frozen protocol](PROTOCOL.md) predates expression aggregation and DE
outcomes. Nine disjoint cohorts contain twelve donors each, six per sham arm,
with the original batch fixed effects. Assignment uses metadata and fixed
hashes, never expression or significance. All released cells use the author's
broad `AIFI_L1 = B cell` annotation, including the five original B/Plasma
subtypes. No subtype is selected after seeing outcomes. Each donor has one
original sample/visit in this release; the original donor, sample-kit, batch,
age-group, sex, CMV and subtype-count metadata remain in `cohorts.json`.

The public CELLxGENE release is marked `is_primary_data=false` because its cells
also occur in other corpus releases. Those overlapping representations are not
additional independent observations. This is our assignment-null experiment,
not a reproduction of the author's aging analysis.

## Implementation repairs

The original CLI rejected this source at its old 1 GiB byte limit. Streamed
source admission now shares the existing 64 GiB fixed-buffer snapshot ceiling.
PCA reference/query snapshot copies use the same ceiling. Count, cell, feature,
aggregate, plan and report limits are unchanged. A regression copies a logical
1 GiB-plus-19-byte source through the actual downstream snapshot paths; another
rejects a sparse logical file above 64 GiB before copying it. Those boundary
files are not HDF5/model qualification fixtures.

Native publication of the complete real H5AD took 302.24 seconds with 817,709,056
peak resident bytes; replay took 309.89 seconds with 881,557,504 bytes on the
physical M4 Pro. The independent SciPy audit compares every raw aggregate,
original axis identity, sample assignment, group membership and cell QC value.
It matches exactly. Neither route constructs a dense cell-by-gene matrix.
Metadata, QC and donor-by-gene aggregates remain resident; this is not a
64 GiB-file, million-cell, GPU or end-to-end performance qualification.

The [separately declared follow-up](FOLLOWUP.md) repairs an unconditional
Cook's-distance requirement. Singleton batches have unit-leverage observations
whose influence diagnostic can be unavailable despite an identified treatment
contrast. Wald/LRT now require that diagnostic only when an influence threshold
is requested. Missing diagnostics remain missing; an explicit unevaluable
influence policy still withholds inference. No normalization, dispersion, count,
support, batch, boundary, test, default or significance threshold was tuned.

The original 27 native runs remain archived. The repaired Wald and LRT methods
run on all nine cohorts as 18 separately labeled `native-fixed` cases. Their
count/design/normalization/trend and final fits are compared exactly against
the baseline; previously available raw probabilities are unchanged. Expanded
test families can change BH values. The source-size fix passed nineteen initial
tests. The count-model repair passed 42 tests in nine suites; the final snapshot
chain passed twelve tests in five suites, including the new byte-boundary test.
These suites overlap and are not an additive count of distinct tests.

## Comparisons and evidence

All 27 original native runs, eighteen repair runs and 27 reference runs
completed. Numerical checks pass for 662,960 retained converged fits and 492,774
native tests across original/follow-up analyses; these are repeated analyses,
not distinct biological genes. The repair restores 85,093 tests per Wald/LRT
method across the nine cohorts. Existing fits and available raw probabilities
are exactly equal in the retained diagnostics, while BH expands to the complete available family.

The following comparison uses the 132,592 gene/cohort hypotheses tested by all
six methods after the repair, with BH recomputed separately per cohort/method:

| Method | BH ≤ 0.05 calls across nine cohorts | Cohorts with any call |
| --- | ---: | ---: |
| Repaired native Wald | 115 | 9/9 |
| Repaired native LRT | 108 | 9/9 |
| Native adjusted QL | 2 | 1/9 |
| edgeR robust QL | 2 | 1/9 |
| limma voom/robust EB | 3 | 1/9 |
| DESeq2 unfiltered Wald | 51 | 9/9 |

These are sham calls. The results expose a Wald/LRT calibration concern;
the small QL call inventory is not sufficient to qualify calibration or choose
a production default. Full reference families also include source genes withheld
by native support rules: DESeq2 has fifty unfiltered calls and 67 default-policy
calls on those original families. Native baseline Wald/LRT had only 47,499 tests
and twenty/sixteen calls because of the diagnostic bug, versus 132,592 tests per
repaired method. The baseline's smaller inventory must not be presented as
better error control.

Each original cohort runs native Gamma-trend NB Wald, LRT and modern adjusted
QL. Pinned references are edgeR 4.10.5 robust QL, limma 3.68.5 voom/robust EB and
DESeq2 1.52.0 Wald, with statmod 1.5.2 and jsonlite 2.0.0. References use exact
independent aggregate counts and the same numeric design. NumPy independently
reconstructs native median-ratio offsets; every native result is checked against
them at a 2e-13 relative/absolute scale tolerance. All original feature rows are
retained before each method's declared count/support/filter policy. DESeq2's
default filtering and explicit unfiltered BH are separate outputs.

[Evidence](evidence/2026-09-10/manifest.json) retains complete original and
repaired probabilities/effects, fitted means and diagnostics, reference tables,
all warnings/failures, full source metadata/provenance, independent aggregate/QC
arrays, requests, logs, runtime identities and both standalone executables.
The external H5AD URL, exact bytes and SHA256 identify the original download.

`families.tsv` and `families.json` report attempted/eligible/tested/withheld
genes, native statuses/boundaries/errors, raw P fractions and BH thresholds.
`joint-families.json` separately recomputes BH on genes tested by all six methods
in each comparison; it does not replace any original method family. Native
checks independently verify source/design identity, NB scores/likelihoods,
conditional information, constrained LR, normal/chi-square/F tails and BH.
QL tail checks condition on its reported denominator DF and verify the family
cap; successful cohort output does not retain its individual moment-quadrature
records, so this is not a new independent moment/DF requalification.

Initial retained execution issues include the old source-size rejection,
missing explicit HDF5 library configuration, a Python archive string-axis loader
error, a runner receipt-key collision and a checker assumption that withheld
fits always expose an effect. Each correction is recorded separately from the
statistical outputs. No failed run is silently replaced in the archive.

## Reproduction

Download the protocol's pinned source and verify its size/SHA256. `prepare.py`
freezes metadata-only cohorts; `audit_counts.py` scans all raw CSR values in
512-cell blocks; `requests.py` writes the 27 original contrasts. Publish and
verify the original source with the native pseudobulk CLI and explicit native
HDF5 library. Keep `report.json` unchanged in that bundle and create a separate
`gzip -n` copy as `source-report.json.gz`, with receipt `native-receipt.json`.

Run `check_source.py --root ROOT`, then `prepare_reference.py --root ROOT`.
Build the inference harness with `build.sh OUT`; use `run_native.py --root ROOT
--binary OUT/independent-null`. The original baseline uses commit
`73aea31c00a7adfc9743b4373efa3f2c17e73b08` for its inference owner. The repaired
owner uses `--output-subdirectory native-fixed` and the separately built binary.
Run `run_reference.py` with explicit `--rscript` and `--r-library` paths.
Drivers checkpoint completed cases and refuse ambiguous unfinished destinations.

After all jobs finish, run `check_models.py` for each `--variant native` and
`--variant native-fixed`, then `check_repair.py`, `summarize.py` and `archive.py`.
All require `--root ROOT`; the archiver also requires a new `--out DIRECTORY`.
`--available` permits explicitly partial inspection during execution.

Nine cohorts from one selected study are not nine independent studies or a
precise universal FDR estimate. Event fractions are descriptive, with no IID
binomial interval. This experiment supplies no alternative-model power,
effect-interval coverage or biological-intervention truth. Production defaults
remain unchanged; broader independent calibration and perturbation prediction
remain part of the full single-cell goal.
