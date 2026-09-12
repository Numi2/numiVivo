# HIRISA integration design: donor, batch and protected preparation

This audit checks the full source-bound PCA mapping against every original
observation's GEO library and preparation metadata: 1,612,594 cells in 131
libraries, five donors, thirteen batch labels and eight condition labels.
It does not run a new correction or claim biological preservation.

The present native integration owner accepts one donor or batch covariate.
Its connectivity admission check protects separation from condition labels;
it does not check arbitrary biological/preparation covariates. This audit is
preparation for multiple-covariate integration, not an implementation of it.

## Identifiability result

Columns below use an intercept and reference-coded categorical indicators.
Rank is checked independently by SVD and column-pivoted, reorthogonalized QR;
weighting libraries by their observed cell counts leaves every rank unchanged.

| Source preparation | Libraries | Cells | Donor + batch rank / columns | Condition dimensions estimable after nuisance |
| --- | ---: | ---: | ---: | ---: |
| Complete atlas | 131 | 1,612,594 | 17 / 17 | 7 / 7 |
| Bcell | 25 | 300,665 | 6 / 10 | 4 / 4 |
| Monocyte | 25 | 319,578 | 5 / 9 | 4 / 4 |
| NK | 25 | 289,717 | 5 / 9 | 4 / 4 |
| PBMC | 26 | 327,062 | 6 / 6 | 2 / 2 |
| Tcell | 30 | 375,572 | 10 / 10 | 4 / 4 |

The complete atlas supports separate donor/batch design columns, but Bcell,
Monocyte and NK subsets do not. Ridge can make their penalized coefficient
solution unique; that does not make the individual nuisance effects identifiable
from the data. Estimable condition main-effect dimensions likewise do not prove
causal identification or protection of other biology.

Across the complete atlas, **Monocyte preparation is exactly equivalent to
membership in these five batches**:

- EXP-00602-PA
- EXP-00633-PA
- EXP-00667-PA
- EXP-00671-PC
- EXP-00679-PC

Every library satisfies that indicator equality. No expression-derived cell
labels were introduced: these are deposited preparation labels, not an assertion
that every selected cell is a biologically pure monocyte.

Condition plus preparation spans ten non-intercept design dimensions on its own,
but only nine additional dimensions after donor/batch nuisance columns. One
protected dimension is already in the nuisance span. The condition/preparation
matrix itself also has one redundant reference-coded column; its baseline rank,
not the raw column count, is used for this comparison.

## Consequence for implementation

Do not add an unconstrained donor-plus-batch correction and claim that the
condition-connectivity test guarantees biological preservation. An explicit
protected-covariate declaration and identifiability report are needed first.
A correction must state how it handles nuisance directions shared with protected
preparation or biology, and retain them or reject the conflicting request rather
than claiming to estimate both independently. The existing donor default and
historical integration outputs remain unchanged.

The next qualified integration change should include real-source crossed and
nested designs, independent numerical reconstruction, preserved original PCA
and RNA counts, and the same sensitive marker/program and rare-cell tests.
This metadata result neither explains all 37/146 historical HIRISA failures nor
turns their failed preservation scores into passes.

## Reproduction and evidence

`audit.py` verifies the original H5AD SHA256 and the PCA plan/receipt hashes from
the frozen full-clustering manifest. It checks every cell's library and source
`geo_enrichment`, rejects inconsistent library preparation labels, and retains
all 131 source sample definitions and cell counts. It then computes every design
rank in all six scopes and verifies the exact preparation/batch identity.

Run with the retained scientific Python environment containing AnnData, h5py and
NumPy. Adapt the explicit paths if restoring elsewhere. The original HIRISA H5AD
and three small provenance files from `numivivo-cluster-export-20260912` are
required; the original data are not copied into Git. `design-audit-initial.json`
retains the first donor/batch/condition-only audit; the final report adds the
protected-preparation checks. Neither phase reads expression outcomes to select
a correction. The archive contains both reports, script and execution logs.
