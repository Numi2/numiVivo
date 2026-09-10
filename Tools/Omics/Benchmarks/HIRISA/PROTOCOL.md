# HIRISA experimental and scale benchmark

Frozen 2026-09-10 at implementation 865f6f6a1963389f4fa4a475fe06b34712123c2c,
before count-value inspection, aggregation, fitting, or prediction. The complete
GEO SOFT metadata, archive listing, author methods, and the HDF5 hierarchy,
shapes and dtypes of the first accession have been inspected. No DE tables,
response signatures or observed count values have been used for selection.

## Source and interoperability

Acquire all 131 individual labeled H5 files from GSE306664, retaining their
original bytes, URLs, file sizes and SHA256 hashes. Keep all source observations
and features, with no subsampling, cell QC selection or expression filtering
at ingestion. The author resource describes over 1.2 million cells; audit the
actual deposited count rather than asserting that scale from the description.
The release contains 105 enriched-population libraries and 26 whole-PBMC
libraries from five donors. It is 10x Flex fixed-RNA, with CellRanger counts
and author cell annotations; it is not an unprocessed sequencing archive.

The first source uses the 10x feature-by-cell CSC layout with uint16 data and
indices, uint32 indptr, and an additional matrix/observations group. Verify
all files, their feature axes, all sparse offsets and indices, integer counts,
cell IDs and annotation lengths. Preserve every observation column and feature
annotation in an AnnData representation or the retained source, documenting
which representation carries each field. Namespacing by GEO accession must
not conceal duplicated biological cells: audit original cell_uuid across all
files and stop for reconciliation if duplicates occur. Never round normalized
values into counts. Never create a dense cells by genes matrix.

Use complete data to exercise native H5AD streaming and independently compare
all pseudobulk counts and per-cell QC against a chunked Python/SciPy calculation.
Report cells, genes, source sparse entries, canonical nonzeros, UMIs, exact
source hashes, elapsed time and peak memory. Admission ceilings alone are not
scale evidence. Retain structural failures and repair the native owner where
required; do not reduce the study to make existing limits pass. Separate
Python format preparation from native ingestion and downstream execution.

## Experimental comparisons and confounding

Keep donor, accession, enrichment population, treatment, experiment batch,
pool, well and chip provenance. Donors are 2616BW, 3283BW, 3491BW, 3955BW and
6811BW. Form all sixteen enriched-population contrasts: Bcell, Monocyte, NK and
Tcell crossed with IFNa, IFNb, IFNg and IFN-L1 versus none. Match controls by
original donor, enrichment population, experiment batch and pool. Require a
unique matched control, never pooling the two distinct Tcell control batches.
GSM9205558 (3491BW Bcell IFNg) is in EXP-00756 without a same-batch enriched Bcell
control: retain its counts and provenance, but withhold that donor pair from
confirmatory DE and prediction. The resulting Bcell IFNg contrast has four
eligible donors; the other fifteen have five, subject to count eligibility.
Batch effects aliased with donor are absorbed by paired donor terms, not
included as redundant columns. Report the exact numeric design and rank.

Use native paired-donor NB Wald, LRT and adjusted QL, median-ratio offsets,
gammaParametric dispersion trend, minimum twenty trend genes, minimum prior
variance 0.25 and two-SD dispersion outlier threshold. Require three paired
donors, ten cells per pseudobulk, total gene count ten and expression in three
pseudobulks. No Cook's exclusions, count replacement, active-donor selection,
effect prior or dispersion fallback. Preserve failures and withheld genes.
Compare pinned edgeR robust QL, limma voom/robust eBayes and DESeq2 Wald using
identical counts, donor design and native offsets. Keep full method-specific
families and a separate jointly tested family, unfiltered BH results, effect
concordance and all warnings. Cross-method agreement is not biological truth
or FDR calibration. Do not promote a production default from this study.

## Prediction and biological contexts

For each enriched contrast, hold out each eligible donor in turn, train on all
remaining eligible donors, and predict the held-out treated expression from
that donor's control. No held-out treated values may select genes, weights,
regularization or methods. Freeze the existing native response model's options
before fitting. Report all folds and compare the identity/no-response baseline
and the training-donor mean response baseline with identical feature families.
If model settings cannot be frozen from existing qualified owners, publish a
separate pre-fit specification before any response fitting.

The PBMC arm has Fresh, culture_IFNa and culture_no_stim libraries. Fresh is not
a 21-hour matched control. Retain all libraries for interoperability/scale;
for a later transfer test use culture_IFNa versus culture_no_stim, matching
donor and original pool and combining technical replicates only within the
same biological condition. Keep author annotation labels as predictions, not
verified cell identities. A cross-preparation transfer analysis must freeze
its label mapping and training/test split before fitting. Never describe
whole-PBMC technical replicate libraries as independent biological donors.

This protocol sets the complete study target. An acquisition audit, one-file
import, or numerical regression is not completion of the native benchmark,
prediction, million-cell out-of-core pipeline or overall single-cell goal.

Sources:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306664
- https://apps.allenimmunology.org/aifi/resources/ifn-response/
- https://apps.allenimmunology.org/aifi/resources/ifn-response/methods/
- https://apps.allenimmunology.org/aifi/resources/ifn-response/analysis/
- https://doi.org/10.64898/2025.12.02.691676
