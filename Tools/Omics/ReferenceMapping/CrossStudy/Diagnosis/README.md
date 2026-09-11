# Why the rare source labels did not transfer

The complete Kang–Ding experiment still **fails both declared transfer targets**.
This post-result diagnosis traces a source-label discrepancy through the original
H5AD and both frozen classifiers. It does not relabel cells, exclude a stratum,
refit a model or revise the gates. All 68,704 query cells remain represented.

## Measured transcript discrepancy

Every one of the 24,673 original Kang barcodes and decoded cell-type labels
matches the prepared input exactly. The literal `Megakaryocytes` label is category
code 7 for all 132 cells. Thus the preparation did not introduce this annotation.

Two explicitly analyst-selected transcript means were scored by the existing
native program owner: HBB/HBA1/HBA2 and PF4/PPBP/ITGA2B/GP9/TUBB1. Each score is
the equal-weight mean of log1p(count / full original RNA library × 10,000).
All members must be present. These post-result gene sets are descriptive RNA
summaries, not independently validated identity, pathway or protein measurements.

| Literal source-labelled class | Cells | HBB/HBA1/HBA2 mean | PF4/PPBP/ITGA2B/GP9/TUBB1 mean | Any of latter five detected |
| --- | ---: | ---: | ---: | ---: |
| Kang `Megakaryocytes` | 132 | 7.635571 | 0.019861 | 5.303% |
| Ding `Megakaryocyte` | 800 | 0.034893 | 2.507374 | 96.875% |

The archive retains every-cell score and count plus summaries for every original
label and library. No cell was selected by a gene-expression threshold. The
large transcript difference raises a concrete source-label concern; matching
names do not establish biological equivalence. These observations do not prove
a replacement identity or the complete biological cause of the discrepancy.

## Frozen classifier trace

Both classifiers recognize their own training source class almost perfectly:
Kang recall is 1.0 and Ding recall is 0.99625. This is an in-sample diagnostic,
not generalization. Both models retained all eight diagnostic genes among their
selected features. A missing panel member does not explain this failure.

For each frozen PCA/scaling/logistic model, the checker reconstructs every query
class logit directly from gene expression and compares every class probability
with the original native report. It then decomposes the mean source-class versus
CD4 T logit margin into all selected gene contributions. CD4 T was chosen after
seeing that it dominates the observed errors; this is an exploratory contrast.

| Frozen model | Training literal-class mean margin | Other-study literal-class mean margin | Difference |
| --- | ---: | ---: | ---: |
| Kang → Ding | 11.451242 | −20.193375 | −31.644618 |
| Ding → Kang | 23.995484 | −5.661905 | −29.657390 |

For Kang → Ding, HBB, HBA2 and HBA1 contribute −9.555738, −5.597256 and −4.406607
to that difference. For Ding → Kang, NRGN, PF4 and PPBP contribute −2.632099,
−2.252231 and −2.088482. Every selected gene, coefficient, mean, signed
contribution and query-cell margin is retained. These are algebraic contributions
to a fixed model contrast, not causal biological effects.

The original six-family macro-F1 scores remain 0.612564 and 0.665731. Recall remains
zero for all 800 Ding and 132 Kang source-labelled megakaryocytes. A label-validity
concern limits interpretation; it does not turn either frozen experiment into a
success. Differences in assay, context, expression representation, taxonomy and
partial label coverage remain relevant to broader transfer reliability.

## Source provenance and literature context

Original Kang input SHA-256:
`e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830`.
The original release is available from [Figshare](https://ndownloader.figshare.com/files/34464122).
Prepared Kang SHA-256:
`1d48c1ff10bcfad7c1bc75863ce702a491c6c71b08e921a5e9754a986a9d19f3`.
The checker independently decodes the legacy compound `obs` and category array;
it checks all category codes are in range and all labels/barcodes agree.

A primary institutional copy of Seninge et al., *VEGA is an interpretable
generative model for inferring biological network activity in single-cell
transcriptomics* (2021; DOI 10.1038/s41467-021-26017-0) has an indexed methods
excerpt: “Megakaryocytes were removed due to uncertainty about their annotation.”
[Institutional article PDF](https://escholarship.org/content/qt0j91d9q3/qt0j91d9q3_noSplash_51497787a7c537c963ad02a20ca2b749.pdf)

Only the search-indexed excerpt was available in this audit. Full PDF retrieval
returned HTTP 429; PMC presented a challenge, and the publisher route failed.
Those routes were not bypassed or repeatedly retried. The retained literature
record states these retrieval limits. The excerpt provides historical context,
not proof of exact source-cell identity or permission to adopt that exclusion.
The local count and original-label checks above are direct, separate evidence.

## Execution and reproduction

Four native score/verify commands completed over both entire sources. Independent
sparse arithmetic agrees with every score to a maximum absolute error of
1.776357e-15; detection counts and full-library totals are exact. Every frozen
query logit reconstructs with maximum absolute error 4.618528e-14. Each full
mean-margin difference equals the sum of gene contributions within 1e-10.

Execution used the already qualified CPU executable from owner revision
`4120f929fcd0e96c8731ee678a320a9b5c7ff11d`, SHA-256
`9d5de14bb4e02f5d6a8db46df95d1bfc9b27205119afca4084139608fb6e0537`.
HDF5 SHA-256 is `a00ffbf8ab94ad81f67231a1ae01df748689e1c35a3615a57f88f4710b8d213e`.
No native source changed or classifier was refitted. The physical Apple M4 laptop
has 24 GiB RAM and macOS 26.6; maximum native RSS was 133,791,744 bytes. This is
an operational CPU measurement, not a speed comparison, GPU or full-app result.
The checker used AnnData 0.13.3.post0, NumPy 2.5.3 and SciPy 1.18.1.

The original prediction freeze remains
`d7cc8f7e72a3918332914256ebeeab207f32443188cdbbfcdd7fa87d99264261`.
Diagnostic native input freeze:
`c1be406588b48c12ec3c92f7ca99573d3c54f8358a246579123041dd1e0eb568`.
The initial informal Python trace failed because `inspect.py` shadowed the
standard library; renaming it to `trace_reference.py` resolved the import. Its
failure and successful trace are retained separately from native checks.

The [protocol](PROTOCOL.md), [native runner](run_programs.py) and
[independent checker](check.py) accept explicit paths with `--help`.
Use the previous [complete cross-study inputs and prediction archives](../README.md#reproduce-and-retain).
The [compact archive](evidence/2026-09-11/manifest.json) retains all native metadata,
plans, receipts and binary outputs, all diagnostic arrays, logs and recipes.
Complete original H5AD bytes remain in the prior input directory and are bound
by size and SHA-256 as explicit replay dependencies; they are not duplicated
inside this diagnostic archive. To restore the two native bundles:

```sh
python3 Tools/Omics/ReferenceMapping/CrossStudy/Diagnosis/archive.py restore \
  --archive Tools/Omics/ReferenceMapping/CrossStudy/Diagnosis/evidence/2026-09-11 \
  --prior-inputs /path/to/numivivo-reference-cross-study-20260911/inputs \
  --out /new/path/to/programs
```

Run the recorded executable's `singlecell-h5ad-programs-verify` on each restored
`kang` and `ding` subdirectory with the recorded HDF5 library configured.
Native replay validates the stored arithmetic and provenance, not cell identities.

## What advances the prediction goal next?

Do not optimize this classifier to force uncertain source labels to agree.
Preserve the two failures and make annotation validity an explicit requirement
for a subsequent frozen transfer experiment. A new benchmark needs documented
compatible label definitions, independently supported rare strata and an untouched
query evaluation, with uncertain and unavailable labels reported explicitly.
Calibration and novel-class rejection still need their own evaluation.

For the strategic perturbation goal, continue the independently collected
GSE181897 input's primary experimental-role resolution before fitting response
models. Its complete count admission already passed; the single-letter B/C
intervention meaning remains unresolved. Neither this transcript diagnosis nor
an external curator's label supplies that missing experimental fact. Reliable
RNA-response transfer, RNA-to-phenotype links and general biological outcomes
remain unqualified; see the [biological prediction assessment](../../../../../Documentation/BiologicalPrediction.md).
