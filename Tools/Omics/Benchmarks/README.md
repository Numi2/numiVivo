# Real experimental single-cell comparisons

The [full-source pancreas benchmark](PANCREAS.md) adds all four Baron human
donors and 14 author-assigned cell types, with native streamed QC/pseudobulk
checks against every original source row and Scanpy. It also retains the rejected
fractional Muraro release and the unresolved streaming-PCA/integration scope.

Completed native comparisons run the actual product commands and replay their receipts.
Python supplies public-data preparation and independent Scanpy/PyDESeq2
references. It is not a dependency of native NumiVivo count execution.

The [complete HIRISA source benchmark](HIRISA/README.md) adds all 131 deposited
libraries: 1,612,594 cells and 3.846 billion count entries, with exact source and
backed AnnData checks. Native release ingestion/reconstruction and independent
verification of every cell/QC value and all library aggregates now pass. Peak
RSS remains 4.43/4.96 GB; metadata/report storage is still resident. Donor/batch
matching is frozen, reference DE runs are complete, and native QL comparisons
remain active. All 79 held-out prediction folds are frozen and unfitted.

## Recorded scope

The [2026-09-09 evidence](evidence/2026-09-09) retains source and executable
SHA256 hashes, exact designs, versions, command timings, peak resident memory and
reference warnings. Downloaded data and large native artifact stores are not
committed. Paths in these measured receipts identify the original local run.

| Public release and declared scope | Result | Remaining boundary |
| --- | --- | --- |
| [Kang 2018](https://doi.org/10.1038/nbt.4042), all 2,651 annotated B cells and all 15,706 source genes | Exact native counts, 16 pseudobulks and Scanpy QC. Eight paired donors, IFNB versus control: baseline 8,894 tested; native NB 5,400 tested and 3,494 support-rank rejections. NB effect correlation 0.99907 versus PyDESeq2. | One cell type and contrast; no calibrated significance or production claim. |
| [PBMC3K](https://scanpy.readthedocs.io/en/latest/generated/scanpy.datasets.pbmc3k.html), all 2,700 cells and 32,738 genes | Annotation preservation, full 2,286,884-nonzero native count/replay route, exact Scanpy QC; normalization max error 8.88e-16. | One library, no donor inference. Legacy source requires current AnnData re-encoding. |
| [Haber 2017](https://doi.org/10.1038/nature24489), all 409 annotated tuft cells and 15,215 source genes | Exact native count/replay and Scanpy QC; normalization max error 8.88e-16. Source dense counts read in 64-row blocks into sparse storage. | Individual donor IDs unavailable in this release; only two batch labels per treatment. No donor-DE claim and no invented donors. |
| [Hagai 2018](https://doi.org/10.1038/s41586-018-0657-2), supplied pertpy release | Ineligible for count DE: fractional X and no raw/counts layer. | Remains ineligible; original author counts are a separate release below. |
| Hagai original E-MTAB-6754, all six mouse unstimulated/LPS6 files | 13,863 cells, 22,048 genes, 32.85M nonzeros; exact source/Scanpy QC and pseudobulks. Native NB/PyDESeq2 comparison on three source-derived donor pairs; 12,426 tested, correlation 0.99919, five expected genes positive in both. | Pairing inference from source prefixes is explicit. Two residual degrees of freedom; substantial p-value differences; six-hour contrast is not the paper's four-hour DE reproduction. |
| [Crowell cortex](Crowell/README.md), all 25,224 deposited post-QC nuclei and 11,076 genes | Eight independent mice; exact counts/QC/offsets, seven native analyses/replays and 42 edgeR/limma/DESeq2 fits. 72,739 tested gene/population pairs. | CPE replication gate retained; eight strict refit-p checks fail with documented stopping-tolerance explanation and no changed BH<0.05 decisions. Strong significance disagreements remain; no FDR calibration claim. |

Dataset download IDs follow the [official pertpy data definitions](https://pertpy.readthedocs.io/en/1.0.2/_modules/pertpy/data/_datasets.html).
Checksums in the scripts bind the exact releases, not every dataset associated
with the corresponding paper. `audit_hagai.py` qualifies only ineligibility of
this supplied file, not absence of count data elsewhere.

## Kang comparison

Both workflows use identical integer pseudobulks and the same count/presence
feature filter, with design `~ donor + condition`. The native method is the
existing median-ratio, moderated log-linear baseline. The reference is PyDESeq2
0.5.4 with a negative-binomial likelihood; it is not R DESeq2, edgeR or limma.
Batch metadata is unreported, so no batch adjustment is asserted.

The declared expected direction check passed all five genes (ISG15, IFIT1,
IFIT3, MX1, OAS1), against a requirement of at least four positive effects in
both models. Effect Spearman correlation was 0.8991, sign agreement 86.42%, and
26 of the top 50 BH-ranked genes overlapped. These are descriptive comparisons;
no parity or competitiveness threshold was set. The substantially different
p-values are retained in per-gene outputs, not used to claim model equivalence.

PyDESeq2's parametric dispersion trend did not converge and it fell back to a
mean trend; this warning is retained in the result. Cook's refitting/filtering
and independent filtering were disabled for this controlled same-filter
comparison. Robust reference sensitivity remains unqualified. Cell counts are
not treated as independent DE replicates.

## Reproduction

Use Python 3.12 and the complete pinned `requirements.txt`. The recorded native
binary was a Swift debug build on Apple silicon, SHA256
`3420335413edc5cbfe3b25c7af7f95731708423dfd0ac9043cc3ba78556e8cc0`.
`omics-source-sha256.json` records the tested owner files. Command timings and
RSS are scope measurements, not optimized throughput or cross-tool performance
comparisons. Reference fit time excludes native import/replay costs and must
not be compared with the native end-to-end command sum.

Download the public source files; the scripts verify SHA256 before analysis:

```sh
curl -fL https://ndownloader.figshare.com/files/34464122 -o kang_2018.h5ad
curl -fL https://ndownloader.figshare.com/files/54169301 -o haber_2017.h5ad
curl -fL https://ndownloader.figshare.com/files/46978846 -o hagai_2018.h5ad
python -m pip install -r Tools/Omics/Benchmarks/requirements.txt
python Tools/Omics/Benchmarks/run_kang.py --binary /path/to/numivivo --source kang_2018.h5ad --out /tmp/kang-result
python Tools/Omics/Benchmarks/run_haber.py --binary /path/to/numivivo --source haber_2017.h5ad --out /tmp/haber-result
python Tools/Omics/Benchmarks/audit_hagai.py --source hagai_2018.h5ad --out /tmp/hagai-audit.json
python Tools/Omics/H5AD/check_public_data.py --binary /path/to/numivivo --full-product --count-analysis --out /tmp/pbmc3k-result
```

Output paths must be new. Run Python without `-O`, since scientific integrity
assertions are enabled by the normal interpreter. Native commands never use
Python for their import, count processing or model fit. Cells-by-genes remains
sparse in both reference and native processing; only the small donor-by-gene
pseudobulk is densified for PyDESeq2.

## Open qualification work

Several independent studies with verified donor structures and raw counts,
additional cell types/perturbations, robust reference sensitivity, explicit
competitive acceptance criteria and significance calibration remain.
[Direct R edgeR/limma-voom/DESeq2 comparisons](../Bioconductor/README.md) now cover
both Kang and Hagai under common count filters and paired designs. Native support
rejections and different multiple-testing families remain explicit.
The [native NB cohort method](../NegativeBinomial/README.md) now has a full Kang
CLI comparison and independent calculation checks. It remains experimental;
old plans retain the log-linear baseline. The count default is
now five million nonzeros; this is still bounded in-memory execution, not
million-cell or out-of-core qualification.

## Original Hagai count recovery

The normalized pertpy release remains ineligible for count DE. Separately,
[the original author release E-MTAB-6754](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6754)
provides integer UMI matrices after the authors' QC/cluster selection. The
[audit](evidence/2026-09-09/hagai-original-audit.json) pins six source files for
mouse sample groups 1–3, unstimulated versus six-hour LPS: 13,863 cells,
22,048 common genes and 32,848,185 nonzeros. This is acquisition and count-axis
evidence. The later NB comparison below records the primary individual table
and deposited sample identities. The streaming route supports this
complete source scope without outcome-based or capacity-driven subsetting.

`audit_hagai_original.py --source-dir /path/to/files --out /tmp/hagai-audit.json`
checks the pinned downloads using one gene row at a time. File URLs and hashes
are in the report. Do not invert the normalized H5AD to reconstruct these counts.

## Full Hagai streaming comparison

`prepare_hagai_stream.py` checks the audited source hashes and converts all six
files, one gene at a time, into a gzip-compressed CSC AnnData file. It also
retains independently accumulated cell totals, detected genes and six sample
pseudobulks. Both sparse index arrays use Int32 because mixed Int32/Int64 arrays
caused the tested SciPy QC path to reject our first prepared file. No cells or
genes are removed; original Ensembl IDs are retained without inferred symbols.

```
python Tools/Omics/Benchmarks/prepare_hagai_stream.py --source-dir /path/to/files --audit Tools/Omics/Benchmarks/evidence/2026-09-09/hagai-original-audit.json --out /tmp/hagai-prepared
numivivo singlecell-h5ad-pseudobulk /tmp/hagai-prepared/prepared.h5ad --plan /tmp/hagai-prepared/plan.json --output /tmp/hagai-native
numivivo singlecell-h5ad-pseudobulk-verify /tmp/hagai-native
python Tools/Omics/Benchmarks/check_hagai_stream.py --prepared /tmp/hagai-prepared --bundle /tmp/hagai-native --out /tmp/hagai-comparison
```

All 13,863 cells, 22,048 features, 32,848,185 nonzeros and 132,271,162 UMIs pass
the source/Scanpy QC and exact six-pseudobulk comparison. Source row membership
and ordered gene identities are also checked. Mitochondrial annotation is
unavailable in this preparation and remains missing. The source uses six-hour
LPS; the original paper's cross-species DE used four hours, so this is not a
reproduction of that DE result. The count-only preparation leaves donor IDs
unset, uses source sample IDs for grouping,
and requests no DE contrast. This qualifies count processing, not independent
replication or expected biological effects.

`check_streamed_de.py` separately compares the streaming Kang NB route with
the previously qualified native resident cohort report. It resolves membership
by sample/barcode/group identity because MEX reorders cells by sample, then
requires every design, numerical and diagnostic field to agree exactly.
This regression preserves the existing Kang evidence and its limitations.

The [2026-09-09 evidence](evidence/2026-09-09-streaming/) records the full
executable SHA, changed source hashes, preparation/source fingerprints,
comparisons, reconstruction checks and unsuccessful preparation attempts.
The full Hagai command took 5.73 seconds with 217,694,208 bytes peak RSS on the
local Mac in one run overlapping other validation work. This includes copying
and hashing the source and serializing the report; it is not a cold-cache or
competitive performance benchmark. The complete executable passed 24 streaming
checks and six existing H5AD CLI checks; 27 import checks, 26 annotation checks
and 18 Swift tests also passed. The existing single-cell CLI suite passed
18 assertions across 24 commands.

## Hagai paired negative-binomial comparison

`run_hagai_nb.py` requires the exact supplementary PDF and SDRF before assigning
donors. Supplementary Table 2 (PDF page 12) explicitly lists three individual
eight-week-old female C57BL/6 mice with single-cell time courses. The SDRF names
the corresponding `mouse1`, `mouse2`, `mouse3` time courses and identifies each
library by an ENA accession. Pairing unstimulated and LPS6 libraries by those
prefixes is an explicit inference from the deposited naming, not a direct
prefix-to-table-row crosswalk. No mapping to table rows 4/5/6 or technical batch
is invented. The [supplement download](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6347972/supplementaryFiles)
and its member PDF are hash-pinned in the design evidence.

Before fitting, the script records a paired `~ donor + condition` design,
median-ratio offsets, parametric NB dispersion trends, no batch adjustment,
and a descriptive acceptance rule: effect correlation at least 0.95, zero native
numerical failures, and at least four of Nfkb2/Nfkbia/Cxcl10/Isg15/Tnf positive
in both methods. All six source samples and all genes remain in native input;
the declared gene-level count/support filter defines inference eligibility.
PyDESeq2 Cooks refitting/filtering and independent filtering are explicitly off.

The [recorded result](evidence/2026-09-09-hagai-nb/report.json) passes those
criteria: 12,673 eligible genes, 12,426 tested and 247 rejected for insufficient
positive-count design support, no numerical failures, effect Spearman
0.99918955, sign agreement 99.46%, and 39 common genes in the top 50 BH lists.
All five preselected response genes are positive in both methods. This is a
second experimental NB comparison after Kang, not repeated-sampling calibration.
PyDESeq2 warns that dispersion estimation with only two residual degrees of
freedom is unreliable. Native and reference p-values differ by many orders of
magnitude; native machine-zero tail probabilities remain visible.
The [dispersion comparison](evidence/2026-09-09-hagai-nb/dispersions.json) finds
a median native/reference final dispersion ratio of 0.37368 across tested genes.
For Nfkb2, final dispersions are 0.00849 versus 0.03021. The different trend and
shrinkage estimates therefore supply a concrete source of uncertainty differences.
PyDESeq2 reports 28 false gene-wise and 12 false MAP convergence flags (zero
false coefficient-convergence flags); these are retained without discarding
their genes or claiming that one estimator is calibrated better.

Independent statsmodels/SciPy checks cover all 12,426 native Wald fits, 32 MAP
objectives, the trend objective, prior variance and BH arithmetic. Maximum
absolute log2-effect error is 1.9911e-5; maximum MAP objective error is 5.12e-11.
These validate calculations conditional on the chosen model, not biological
truth or FDR calibration.

```
python Tools/Omics/Benchmarks/run_hagai_nb.py --binary /path/to/numivivo --prepared /tmp/hagai-prepared --sdrf /path/to/E-MTAB-6754.sdrf.txt --supplement /path/to/NIHMS79113-supplement-Supplementary_Information.pdf --out /tmp/hagai-nb
python Tools/Omics/NegativeBinomial/check_cohort.py --cohort-report /tmp/hagai-nb/native/report.json --counts /tmp/hagai-nb/reference-pseudobulk.tsv --out /tmp/hagai-independent
```
