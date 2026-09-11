# Explicit feature panels and complete Kang–Ding transfer

Native reference mapping now accepts different complete RNA gene universes when
an explicit training feature panel is declared. The complete Kang–Ding experiment
passes independent numerical checks for all 68,704 query cells, but **fails both
declared coarse-family transfer targets**. Good aggregate accuracy hides complete
megakaryocyte recall failure in both directions. This does not qualify reliable
cross-study annotation or authoritative biological identities.

## Native input contract

Add `featurePanel` inside `reduction.pca` of a reference-fit plan. It contains
unique exact feature IDs, not display-name aliases or genes inferred from query
expression. Every listed gene must exist in training and query. The gene order in
the list does not affect selection; original source order is retained.

```json
"reduction": {
  "normalizationTarget": 10000,
  "pca": {
    "components": 2,
    "highlyVariableFeatures": 3,
    "featurePanel": ["CD3D", "MS4A1", "NKG7"]
  }
}
```

This fragment illustrates placement; a valid panel must contain at least the
requested number of components and enough genes with defined dispersion. The
complete reproducible plans in the evidence contain all 14,976 declared IDs.

All source genes contribute to each cell's RNA-library total before log1p
normalization. HVG mean bins, dispersion normalization and ranking use only panel
genes. Counts, feature identities and raw moments remain available for all source
genes; genes outside the panel have `meanBin: -1`, null normalized dispersion and
`selected: false`. PCA and classifier parameters use training cells alone.

Reference mapping requires every panel gene, including genes not selected as HVGs.
Other query genes may be absent, different or additional; all measured query counts
contribute to that query cell's denominator. No missing panel gene is filled with
zero. Organism, namespace, count units, cell-overlap and query-label safeguards
remain required. Differing measured universes can still create normalization and
assay differences; the declared panel does not establish biological comparability.

Panel selection is shared by resident and streamed sparse reduction. The different
gene-universe admission described here belongs to `singlecell-reference-map`.
The standalone file-PCA query workflow retains its own input contract.

The optional field is absent in historical plans and outputs. Omitting it retains
the complete-reference-gene-universe requirement and historical numerical bytes.
Uniform kNN remains the default classifier; the experiment explicitly selects
[balanced multinomial logistic regression](../Logistic/README.md). Existing query
score, projection and classifier work limits apply. Panels permit at most 200,000
IDs and at least the requested component count, subject also to serialized-plan
limits. Reference fit CLI plans now use the same 2 MiB bound as reference artifacts
and native file-PCA plans.

## Complete sources and frozen evaluation

The [protocol](PROTOCOL.md) preceded these fits and scores. Both cohorts had
already been inspected for integration, and Baron informed classifier development.
This is a development transfer assessment between separately collected studies.

The existing [Kang source](../../Reduction/FILE_INTEGRATION.md) supplies all 24,673
cells, 15,706 symbols, eight donors and both control/stimulated conditions.
The existing [Ding source](../../Reduction/DING_INTEGRATION.md) supplies all 44,031
UMI-method cells, 33,694 features and twelve libraries in two biological samples.
Ding donor identities are unreported. Its read-count Smart-seq2 cells remain in
the original source provenance and outside this previously declared UMI input.

Of the Kang symbols, 14,976 occur exactly once in Ding's source symbol column;
19 are ambiguous and 711 are absent. Only unique exact matches enter the panel.
Ding receives a separate native feature-ID column for these matches; original
feature IDs and a complete reversible table are retained. No synonym resolution,
duplicate merging or expression-based feature matching is performed. Every count
and both original axes survive preparation exactly.

| Input | Cells | Original genes | Nonzero counts | Total UMIs |
| --- | ---: | ---: | ---: | ---: |
| Kang, training and query | 24,673 | 15,706 | 14,184,532 | 38,667,346 |
| Ding, complete query | 44,031 | 33,694 | 38,102,066 | 92,558,353 |
| Ding, all assigned training labels | 29,411 | 33,694 | 27,449,248 | 68,718,701 |

Ding's 46 Unassigned and 14,574 unavailable labels cannot supervise a fit, but all
of these cells receive retained query predictions. Scored label coverage is
29,411/44,031 (66.80%). The original eight/nine label vocabularies remain unchanged
in native models. Cytotoxic T cell and CD8 T cells are not treated as identical
fine labels; pDC has no separate Kang label. Every original rectangular confusion
table is retained, including unavailable source labels.

Only evaluation groups source subclasses into six declared broad families: T, B,
NK, monocyte, dendritic and megakaryocyte. This coarsening does not establish fine
identity agreement or unknown-class detection. The descriptive target required
balanced accuracy and macro-F1 at least 0.8 in both directions, with each supported
family recall at least 0.5. Cells and technical libraries are not independent donors.

## Results and failures

| Transfer | Query cells | Assigned labels | Family accuracy | Balanced accuracy | Macro-F1 | Training-majority accuracy | Declared target |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Kang → Ding | 44,031 | 29,411 | 0.818945 | 0.654679 | 0.612564 | 0.534494 | FAIL |
| Ding → Kang | 24,673 | 24,673 | 0.900742 | 0.716961 | 0.665731 | 0.521177 | FAIL |

Both transfers have **zero megakaryocyte recall**: all 800 Ding and all 132 Kang
source-labelled cells are missed. Kang→Ding dendritic recall is 0.4860 and NK
precision is 0.4193. Ding→Kang dendritic precision is 0.2853 despite recall 0.9546.
These asymmetric false positives and misses remain visible alongside aggregate
accuracy. Per-family precision/recall and all 28 original query-library summaries
are retained. No best library, condition or donor subset replaces the full result.

All eight native fit/map/reconstruction commands pass. Scanpy HVG selections are
exact in both directions. Maximum independent PCA covariance residual is
1.18145e-11; every query projection agrees within 4.45e-14 using the full original
library totals. Explicit softmax agrees within 1e-12. Tighter independent
scikit-learn fits differ in probability by at most 1.88893e-5, below the fixed
1e-4 bound, with no predicted-label disagreements. Scale, class weights, objective
and gradient checks pass. Repeated scoring reproduces all five output files
byte-for-byte.

Eighteen Swift tests in five suites and seventeen native lifecycle commands pass.
They cover panel-order invariance, full-library normalization, missing unselected
panel genes, additional/absent nonpanel genes, query feature reordering, empty
libraries, prohibited labels, namespaces, work limits and oversized plans.
Default kNN model/report bytes match the previous qualified native executable.

Retained failures include the initial 128 KiB CLI rejection of the 277,706-byte
real plan, a checker environment without Scanpy, and an AnnData fixture assignment
that used an unsupported intermediate LIL matrix. The fixture now converts its
private intermediate to CSR before assignment. Numerical/biological thresholds,
source data and frozen fit plans were not changed in these repairs. Scikit-learn
warnings for library strata missing a predicted family are retained; family
macro-F1 always uses all six declared families.

The result identifies a real transfer limitation. It does not identify its cause
or validate a repair. Assay/context differences, limited source-label coverage,
training representation and label definitions require separate investigation.
Calibration, reliable rare-class mapping, prospective study validation and
biological identity evidence remain open.

The subsequent [post-result diagnosis](Diagnosis/README.md) confirms exact original
Kang label preservation and a sharp measured RNA discrepancy between the literal
rare classes. Both frozen classifiers reconstruct exactly in gene space; neither
was refitted. This identifies a source-label concern and its algebraic contribution
to the failed margins, without establishing replacement identities or changing
either original failure.

## Reproduce and retain

`prepare.py` accepts the existing complete source-linked Kang and Ding prepared
directories and freezes count-preserving inputs plus exact feature identities.
`run.py` fits both directions and freezes verified native bundles before scoring.
`score.py` performs sparse independent checks and then the declared diagnostics.
`regression.py` compares native lifecycle behavior and historical default bytes.
Use each script's `--help` for required paths; no script downloads or selects a
smaller query cohort. Build with `Tools/Omics/H5AD/build.sh OUT --with-cli`, then
run `bash Tools/Omics/ReferenceMapping/CrossStudy/test.sh OUT` with the qualified
`NUMIVIVO_HDF5_LIBRARY` set.

The actual executable SHA-256 is
`9d5de14bb4e02f5d6a8db46df95d1bfc9b27205119afca4084139608fb6e0537`.
All 100 compiled source identities match the published owners. Build/tests used
the physical M4 Pro Mac mini; full transfers used the laptop's physical Apple M4,
24 GiB memory and macOS 26.6 because remote disk capacity was insufficient.
Maximum native RSS was 510,066,688 bytes. This is a CPU operational measurement,
not comparative speed, GPU or full-application qualification. Independent tools
were Scanpy 1.12.4, AnnData 0.13.3.post0, NumPy 2.5.3, SciPy 1.18.1 and
scikit-learn 1.9.0. HDF5 SHA-256 is
`a00ffbf8ab94ad81f67231a1ae01df748689e1c35a3615a57f88f4710b8d213e`.

The input freeze is
`0e6ea669286b6eef3ffbbe7e39e2223f2e5fa26e8a5b8c41a301f289692f70db`;
the prediction freeze is
`d7cc8f7e72a3918332914256ebeeab207f32443188cdbbfcdd7fa87d99264261`.
[Compact evidence](evidence/2026-09-11/manifest.json) retains complete per-cell
probability comparisons, source/class confusion tables, all library metrics,
plans, source identities, logs, failures, tests and executed recipes. Full native
bundles remain locally under
`/Users/home/numivivo-reference-cross-study-20260911/native-qualified`: two verified
deduplicated tar archives preserve 633,223,424 logical bytes in 284,884,173 stored
bytes. Restore with the shared [bundle utility](../Logistic/bundles.py), using the
recorded executable for native verification. The full bundles were not copied to
the nearly full Mac mini; compact evidence and code are published there.
