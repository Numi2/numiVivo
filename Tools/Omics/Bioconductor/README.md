# Direct Bioconductor pseudobulk comparisons

This benchmark runs native NumiVivo NB results against R edgeR QL, limma-voom
and DESeq2 on two existing experimental contrasts. It uses all source cells and
genes before the predeclared count filter. Only donor-level pseudobulk is dense.

- Kang: all 2,651 B cells, 15,706 source genes, eight paired human donors,
  IFNB stimulation versus control; 8,894 eligible genes.
- Hagai: all 13,863 deposited mouse cells, 22,048 source genes, three paired
  individuals, six-hour LPS versus unstimulated; 12,673 eligible genes.
  The deposited sample-prefix pairing inference and primary individual-table
  provenance remain as documented in the [source benchmark](../Benchmarks/README.md).

`prepare.py` independently reaggregates archived H5AD counts and checks every
pseudobulk entry, including filtered-out genes, against the native report. It
also checks the earlier independent eligible-gene reference table. It exports
the exact native observation order, design matrix, contrast, library totals and
size factors. Donor fixed effects and unreported batch metadata are preserved.
It does not infer batches or substitute cells for biological replicates.

## Reference policies

The native reports use the explicit `gammaParametric` dispersion trend; this
does not change NumiVivo's default trend. All methods receive the same count
prefilter: total count at least ten and expression in at least three pseudobulks.
The exact native numeric design matrix supplies the paired contrast.

Two modes separate normalization from other modeling differences:

1. `native_size_factors`: all methods use the native median-ratio factors.
   edgeR/voom effective library sizes multiply those relative factors by the
   geometric mean raw library size. This common scale changes only the
   intercept and keeps voom's library sizes in count units.
2. `package_normalization`: edgeR/voom use TMM; DESeq2 estimates its own size
   factors on the same eligible genes. Full raw library sizes remain supplied
   to edgeR/voom even after filtering.

edgeR uses robust dispersion estimation and robust quasi-likelihood F tests.
limma uses voom precision weights and robust empirical-Bayes moderated t tests.
DESeq2 uses parametric dispersion fitting, Wald tests, no coefficient prior and
no automatic count replacement. Its main comparison disables independent
filtering and Cook's suppression; the ordinary `results()` policy (default
alpha 0.1) is exported separately. No fitted log-fold-change shrinkage is used.
These likelihoods, moderation methods and test distributions remain distinct.

Official references: [edgeR manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf),
[DESeq2 workflow](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html),
[limma](https://bioconductor.org/packages/limma/).

## Results and the coverage gap

All twelve reference fits complete without model warnings. DESeq2 reports no
coefficient nonconvergence. Current native publication/reconstruction succeeds
and reproduces both previous gamma-trend reports exactly. Independent BH checks
pass on the reported p-values. Detailed tables and diagnostics are retained in
[evidence](evidence/2026-09-09/source-state.json).

The main limitation is native coverage: the positive-count support rank gate
withholds inference for **3,494 Kang genes** and **247 Hagai genes**. R returns
results for all eligible genes. Thus effect correlations below cover only
5,400 Kang genes and 12,426 Hagai genes. They do not qualify rejected genes.

| Fixed native size factors | Kang effect Spearman | Hagai effect Spearman |
| --- | ---: | ---: |
| edgeR QL | 0.99017 | 0.99471 |
| limma-voom | 0.98971 | 0.99470 |
| DESeq2 | 0.99853 | 0.99981 |

The ordinary native BH family excludes support-rejected genes, whereas R's main
family contains all eligible genes. Raw significance counts therefore use
different multiplicity families. `compare.py` retains those counts and also
recomputes BH separately on the same jointly testable intersection:

| Fixed factors, joint-family BH < 0.05 | Kang (5,400 genes) | Hagai (12,426 genes) |
| --- | ---: | ---: |
| NumiVivo | 850 | 4,827 |
| edgeR QL | 1,054 | 4,955 |
| limma-voom | 1,009 | 5,043 |
| DESeq2 | 700 | 4,487 |

This intersection calculation is a descriptive sensitivity check, not a change
to the published native inference. With fixed native factors, R calls 35–160
Kang genes and 40–68 Hagai genes significant among the native-unavailable genes,
depending on method. Those are reference calls, not verified biological truths.
All five previously declared response genes per study have positive effects in
all references. Large per-gene disagreements and extreme/zero p-values remain
visible. Package-normalized modes, top-100 overlaps, default DESeq2 filtering,
convergence and dispersion-outlier counts are recorded in the comparison JSONs.

## Reproduce

The checked environment is R 4.6.1 / Bioconductor 3.23, edgeR 4.10.5, limma
3.68.5, DESeq2 1.52.0, statmod 1.5.2 and jsonlite 2.0.0. `run.R` rejects different
versions of these five packages. The complete installed package lock and
`sessionInfo()` records are included. The isolated library used here is
`/Users/home/numivivo-r-library-20260909`; it is reference tooling, not a native
NumiVivo runtime dependency. Requalify package changes explicitly.

```sh
python prepare.py --report /native/report.json --counts /reference/reference-pseudobulk.tsv --expected GENE1 GENE2 --out /new/input
R_LIBS_USER=/path/to/qualified-library Rscript run.R /new/input /new/reference
python compare.py --input /new/input --reference /new/reference --out /new/comparison.json
```

The native report must be in its archived H5AD bundle with `original.h5ad` and
`plan.json`. Native input, reference input and generated tables are hash-bound
in the evidence. Runtime measurements have different scopes: R timings cover
fits on prepared pseudobulks, while native publication includes HDF5 scans and
report creation. They are not end-to-end performance comparisons.

Initial preparation errors (eligible-only tables treated as full axes, then a
NumPy bound-property typo) are retained. The first BH verification exposed
precision loss in pandas' default text float parser; `float_precision=round_trip`
fixed it without relaxing tolerances. Old native receipts correctly rejected the
changed implementation; fresh bundles were published and verified. No stale
receipt was rewritten to manufacture a pass.

## Remaining work

The next native inference investigation is the support-rank coverage gap:
distinguish estimable treatment contrasts with boundary nuisance coefficients
from genuinely unidentifiable or separated effects. Any expanded fit needs
explicit dispersion, degrees-of-freedom, uncertainty and calibration evidence;
removing the existing rank gate alone is not a repair. More independent studies,
cell types, null/calibration experiments, robust nuisance handling and prospective
prediction remain. These two study comparisons do not establish calibrated FDR,
causal validity, general competitiveness or production qualification.

The independent `audit_support.py` algebra check further separates the gap:
3,429 rejected Kang genes have full positive-count design rank after excluding
zero-total donor pairs while retaining at least three active donors. This
identifies candidates for a properly derived boundary-nuisance treatment; no
reduced fit or inference is qualified by that rank check. For Hagai, 169 rejected
genes would retain only two donors, below the declared three-donor minimum,
and 78 remain rank deficient with all three donors. Neither category authorizes
relaxing the inference or replication gates. Per-gene classifications are retained.

Subsequent native work implements an explicit experimental
[active-donor profile policy](../NegativeBinomial/ActiveDonor/README.md), with
independent fits and expanded-family BH checks for all 3,429 Kang candidates.
That policy preserves the original default reports and these original R fits.
Its new descriptive comparisons account for the changed native tested family;
they do not establish selection-adjusted FDR or equivalence to full-donor R methods.
