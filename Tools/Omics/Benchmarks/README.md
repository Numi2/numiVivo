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
Native negative-binomial DE is the next algorithmic step. The count default is
now five million nonzeros; this is still bounded in-memory execution, not
million-cell or out-of-core qualification.
