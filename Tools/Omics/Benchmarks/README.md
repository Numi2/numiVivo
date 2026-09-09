# Real experimental single-cell comparisons

These scripts run the actual native product commands and replay their receipts.
Python supplies public-data preparation and independent Scanpy/PyDESeq2
references. It is not a dependency of native NumiVivo count execution.

## Recorded scope

The [2026-09-09 evidence](evidence/2026-09-09) retains source and executable
SHA256 hashes, exact designs, versions, command timings, peak resident memory and
reference warnings. Downloaded data and large native artifact stores are not
committed. Paths in these measured receipts identify the original local run.

| Public release and declared scope | Result | Remaining boundary |
| --- | --- | --- |
| [Kang 2018](https://doi.org/10.1038/nbt.4042), all 2,651 annotated B cells and all 15,706 source genes | Exact native counts, 16 pseudobulks and Scanpy QC; normalization max error 8.88e-16. Eight paired donors, IFNB versus control, 8,894 eligible/tested genes. | One cell type and contrast; no multi-study competitiveness or native NB claim. |
| [PBMC3K](https://scanpy.readthedocs.io/en/latest/generated/scanpy.datasets.pbmc3k.html), all 2,700 cells and 32,738 genes | Annotation preservation, full 2,286,884-nonzero native count/replay route, exact Scanpy QC; normalization max error 8.88e-16. | One library, no donor inference. Legacy source requires current AnnData re-encoding. |
| [Haber 2017](https://doi.org/10.1038/nature24489), all 409 annotated tuft cells and 15,215 source genes | Exact native count/replay and Scanpy QC; normalization max error 8.88e-16. Source dense counts read in 64-row blocks into sparse storage. | Individual donor IDs unavailable in this release; only two batch labels per treatment. No donor-DE claim and no invented donors. |
| [Hagai 2018](https://doi.org/10.1038/s41586-018-0657-2), supplied pertpy release | Ineligible for count DE: fractional X and no raw/counts layer. | Original counts must be acquired. No rounding or inverse-normalization reconstruction. |

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
competitive acceptance criteria and R edgeR/limma/DESeq2 comparisons remain.
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
evidence, not a completed donor-DE benchmark. Donor mapping must be finalized
from original sample metadata. The streaming route below now supports this
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
reproduction of that DE result. Individual donor mapping requires Supplementary
Table 2, identified by the [primary methods](https://www.hagailab.org/wp-content/uploads/2019/06/Hagai-Nature-2018.pdf).
The preparation leaves donor IDs unset, uses source sample IDs for grouping,
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
