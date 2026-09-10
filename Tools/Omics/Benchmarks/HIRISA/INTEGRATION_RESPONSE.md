# Complete-cohort integration response diagnostic

The subsequent [full-cohort program diagnostic](INTEGRATION_PROGRAMS.md) fails
three control-sensitive comparisons for native integration and all three Harmony
seeds, with 18/32 comparisons insufficiently sensitive. Earlier coarse results
remain valid at their stated scope; complete biological preservation is unqualified.

The original **1,612,594 cells in 131 libraries** have a measured unintegrated
baseline, an exact identity control, and a library-mean-erasure control. The
complete native seed-7 publication now passes all frozen mean-response and
classification preservation margins; native replay and method comparison have
separate status records.
This evaluator adds a bounded experimental-response diagnostic; it does not
complete biological preservation or replace the previous Kang NK-cell failures.

## Frozen comparison

`check_integration_response.py` checks every source barcode, original accession,
donor, treatment, enrichment and batch/pool against the frozen native metadata
and original `design.json`. It reads every 20-PC score record in at most
8,192-row batches, checking coordinates and finite values. Only one uint16
library code per cell persists, plus library means and 20-by-20 second moments.
There is no resident full score matrix or cells-by-genes matrix.

The original sixteen contrasts yield the same **79 donor-held-out folds** as
the matched DE/prediction design. **1,275,710 distinct cells** enter these
classifiers. The **336,884 other cells** remain in the complete moments and
conditional donor-scatter summaries; they are not silently recast as matched
treatment examples. The unmatched Bcell IFNg library and whole-PBMC preparation
differences retain their original exclusions. Enrichment is an experimental
preparation label, not an assertion that every cell has that identity. Author
cell-type predictions are not classifier targets or verified ground truth.

Each linear classifier minimizes equal-library-weighted mean squared error
with ridge one, targets minus/plus one, an unpenalized intercept and threshold
zero. Centering and scaling use training libraries only. Class-balanced held-out
accuracy gives each control/treated class equal weight, then gives donors equal
weight within each contrast. The evaluator uses small sufficient-statistic
normal equations; it does not fit a giant dense training matrix.

Each matched treated-minus-control PC mean response is also compared with its
unaltered baseline, reporting cosine, projected gain and relative Euclidean
drift. A baseline norm at or below 1e-10 is unavailable. Frozen engineering
margins allow mean accuracy loss at most 0.05, individual fold loss at most 0.10,
and donor-mean relative response drift at most 0.25. These are not biologically
calibrated thresholds, confidence intervals, hypothesis tests or FDR control.
Conditional donor scatter uses all libraries within enrichment/treatment strata
and equal donor/library weighting; confounded donor/batch variation prevents
interpreting its reduction as pure technical-batch removal.

The protocol was frozen before baseline and integrated-response inspection.
Baseline admission was subsequently tightened to require the original raw-PCA
score hash and all source bindings, rejecting a corrected result reused as a
baseline. The earlier baseline remains archived. Re-execution under the final
code produced exactly the same moments, classifiers and descriptive results;
no arithmetic, margins or design choices changed.

## Complete baseline and controls

| Enrichment | Treatment | Baseline balanced accuracy | Erasure loss detected (>0.05) |
| --- | --- | ---: | --- |
| Bcell | IFNa | 0.8704 | yes |
| Bcell | IFNb | 0.8790 | yes |
| Bcell | IFNg | 0.9448 | yes |
| Bcell | IFN-L1 | 0.6524 | yes |
| Monocyte | IFNa | 0.8642 | yes |
| Monocyte | IFNb | 0.8858 | yes |
| Monocyte | IFNg | 0.9847 | yes |
| Monocyte | IFN-L1 | 0.5291 | insufficient |
| NK | IFNa | 0.9135 | yes |
| NK | IFNb | 0.9322 | yes |
| NK | IFNg | 0.5196 | insufficient |
| NK | IFN-L1 | 0.5312 | insufficient |
| Tcell | IFNa | 0.8219 | yes |
| Tcell | IFNb | 0.8477 | yes |
| Tcell | IFNg | 0.6272 | yes |
| Tcell | IFN-L1 | 0.5342 | insufficient |

The complete native result has maximum contrast-mean balanced-accuracy loss
0.002546 and maximum individual-fold loss 0.019370. Mean paired response drift
ranges from 0.004105 to 0.138859 across the sixteen contrasts, within the frozen
0.25 margin. Conditional donor-associated variance fractions decrease in all
23 measured enrichment/treatment strata (unweighted mean 0.079429 to 0.020288).
This descriptive reduction includes confounded donor/batch variation. These
coarse diagnostics do not qualify marker programs or rare-cell identity. The
initial remote evaluation failed before metrics because the frozen design file
was absent; copying that exact file and verifying its hash resolved the failure
without changing the evaluator, protocol, margins, baseline or native outputs.
Both attempts are preserved in the full-result report archive.

The identity control reproduces all fold results exactly, with zero response
drift. Subtracting each original library mean erases every paired mean response
(relative drift one in all sixteen contrasts) and returns balanced accuracy
0.5 in every fold. The control analytically sets centered library means to zero
and adjusts second moments; it is a deliberate adversarial control, not an
integration method. Four baseline contrasts are too close to chance for the
frozen classification loss test. Passing a preservation margin there
does not establish classification sensitivity. The mean-response control is
sensitive in all sixteen contrasts, but only tests this coarse response moment.

The final baseline and both controls each scan the complete cohort in roughly
nine seconds locally; maximum measured RSS across them is **102,596,608 bytes**.
This includes source/score hashing and identity checks, but excludes PCA,
integration and count processing. It is not a native/scverse performance claim.

Eight focused tests pass with NumPy 2.5.3 and scikit-learn 1.9.0. Dense weighted
`StandardScaler`/`Ridge` fits independently reproduce every synthetic fold's
coefficients and confusion counts across shuffled library/tile boundaries.
Other checks cover dense moments/scatter, both controls, negligible-response
unavailability, malformed matrix coordinates/nonfinite values/truncation,
source identity mismatches, invalid matched designs and baseline binding errors.
Synthetic tests qualify evaluator arithmetic, separately from the complete
experimental baseline above. No production Swift code changed.

## Evidence and reproduction

The [archive](evidence/2026-09-10-integration-response) retains the exact protocol,
both evaluator versions, original and final baselines, all 79 fold coefficients
and counts, moments for every library, all controls, tests, versions and the
frozen native evaluation coordinator. Its manifest SHA256 is
`aee509905f480ed90745e7ab917a134713a75d37421bdab5de5732b68f3e33c2`.
The archived coordinator waits for the full native publication, validates its
receipt/score/metadata and original input-PC identities, then evaluates that
output with the same frozen baseline and protocol. It does not launch a second
integration run. Native replay and the streamed independent numerical checker
remain separately recorded requirements.

Use the pinned HIRISA requirements; the independent tests additionally require
`scikit-learn==1.9.0`. Set BLAS thread counts to one as in archived commands.

```sh
python test_integration_response.py
python check_integration_response.py \
  --source /study/hirisa.h5ad --metadata /study/pca/metadata.json \
  --design /study/design.json --protocol /study/response/protocol.json \
  --scores /study/pca/scores.bin --out /study/response/baseline.json
python check_integration_response.py \
  --source /study/hirisa.h5ad --metadata /study/native/metadata.json \
  --design /study/design.json --protocol /study/response/protocol.json \
  --scores /study/native/scores.bin --baseline /study/response/baseline.json \
  --out /study/response/native.json
```

Full-cohort transductive PCA/integration precedes these held-out classifiers.
These results do not qualify prospective prediction, marker/program gradients,
rare-cell preservation, learned reference annotation, or unseen perturbations.
A full independent integration-method comparison and the broader twelve-item
goal remain incomplete.
