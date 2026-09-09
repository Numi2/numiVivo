# Frozen PCA query qualification, 2026-09-09

All four predeclared Baron donor-held-out folds are used. Every one of the
8,569 original human cells appears in exactly one query fold, with all 20,125
source genes retained for library normalization. Each training fold selects its
own 2,000 HVGs and fits 20 PCs; no query fitting or query label mapping occurs.

Native commands, elapsed seconds, maximum resident bytes and exit codes are in
`baron/commands.json`. Each fold fits standalone PCA, publishes and reconstructs
the query bundle, then fits/maps the existing label-reference route. All 20
commands pass. Native execution used the Mac mini Apple CPU and isolated HDF5
runtime recorded in `host.txt`; external Scanpy/SciPy checks ran on the laptop.

`baron/` archives every query artifact and complete reference PCA state except
source H5AD files. Decompress `.gz` files to their original names before running
the reference tools. Original query/training H5AD files remain at the paths in
`checks.json`, whose hashes were checked on both hosts. They derive from the
public Baron GSE84133 source and predeclared native projections described in
`../../../ReferenceMapping/README.md`. No measured counts were synthesized.

Every projected score agrees with independent Scanpy normalization and SciPy
sparse matrix multiplication within 2.85e-14 absolute error. All binary coordinates,
source cell/feature IDs, QC totals, selected entry counts and work counts are
checked. Every training score/loading/HVG/statistical result matches the previously
published reference fit. The full legacy query reports, including neighbors,
votes and candidate labels, remain byte-for-byte unchanged. Baseline files were
compared to their committed compressed evidence; decompressed hashes are recorded.
The Python comparison is reproducible with `check_pca_query_baron.py` and the
input/prior paths in its invocation and `checks.json`.

`fixtures-complete/` retains 22 command results including 13 expected rejections:
namespace (including absent training namespace), count units, organism, query
labels, work budget, missing genes, overlapping cells, no-overwrite, and rehashed
score coordinates/values, report counts and training loadings. CSR/CSC/reordered
features and empty libraries are compared independently. `fixture-inputs/` allows
repeating this gate; `legacy-pca/` retains the existing standalone PCA regression.
`tests.log.gz` records 71 tests in 20 suites, including complete score-record checks
across 16 MiB mapping boundaries and rewinds. `release.log.gz` preserves the build
and warnings. No failing native build/test was observed in this qualification.

Query publication peaks at 157-170 MB, including reference reconstruction; these
are end-to-end lifecycle measurements, not isolated kernel timings. Query scratch
is mapped in at most 16 MiB windows. Resident metadata/QC and the training fit still
contribute to peak memory. This benchmark does not establish million-cell scaling,
parallel/Metal performance, corrected donor biology or perturbation prediction.
