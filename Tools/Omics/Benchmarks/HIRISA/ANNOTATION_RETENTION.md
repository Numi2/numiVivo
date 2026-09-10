# Full-cohort annotation retention

**FAIL on the declared complete retention gate for native seed 7 and all three
Harmony references.** Native integration exceeds the loss margins in 37 of 146
supported, control-sensitive label/stratum comparisons, including ten of 29 rare
comparisons. These are losses of recoverable author annotations under the fixed
decoder; the labels are not biological ground truth and the result does not prove
causal cell-type loss.

The [protocol](ANNOTATION_RETENTION_PROTOCOL.md) was frozen after the earlier
HIRISA program, response, clustering and preparation-transfer outcomes were
known, before this annotation diagnostic. No model, label inclusion, margin or
integration setting was selected from these results. Earlier program-preservation
passes and failures remain valid at their stated scopes; they do not substitute
for this separate annotation endpoint. No integration model was refitted.

## Complete cohort and controls

All **1,612,594 original cells**, 131 libraries, five donors, 23 original
preparation/treatment strata and 31 exact author `celltype.l2` labels are retained.
Every cell is a query exactly once across 114 donor-held-out folds; one of the
possible 115 stratum/donor combinations has no source cells. Technical libraries
are pooled within donor, never promoted to independent biological replicates.
No confidence filter, RNA threshold or annotation relabelling is applied.
`Doublet` and very small label groups remain in the ledger, not silently removed.

There are 471 observed label/stratum comparisons. Of these, 330 are operationally
rare: below 1% of source cells within that stratum, determined from metadata only.
Only 163 comparisons have at least twenty training and twenty query cells in
every positive-query donor fold. Seventeen of those lack the declared baseline
advantage over the erasure control, leaving 146 supported and sensitive
comparisons, including 29 rare ones. **308 comparisons lack sufficient support;
17 more have sufficient support but insufficient control sensitivity.** All are
measured and retained; none qualify by averaging only successful donor folds.
Even identity cannot qualify the complete gate with these missing sensitivities.

The fixed class-balanced, one-hot ridge decoder uses all twenty PCs, training-only
weighted scaling, ridge one and an unpenalized intercept. Only other donors enter
training. Classes absent from training retain zero query recall. Exact maximal
score ties split votes uniformly, avoiding an arbitrary class-order advantage.
The all-zero-PC control yields uniform predictions over training-supported
classes and fails all 146 sensitive comparisons. Exact identity reproduces every
fit, confusion matrix and metric. These diagnostic scores are not calibrated
probabilities or a native annotation product.

## Complete comparison

The primary limits are mean donor recall loss at most 0.05 and loss at most 0.10
in every donor fold, conditional on sufficient support and control sensitivity.
These engineering margins do not establish biological utility. The
[complete comparison TSV](ANNOTATION_RETENTION_RESULTS.tsv) retains all 1,884
candidate/label/stratum rows, including insufficient and improved outcomes.

| Candidate | Supported, sensitive comparisons | Failures | Rare supported, sensitive | Rare failures | Complete gate |
| --- | ---: | ---: | ---: | ---: | --- |
| native | 146 | 37 | 29 | 10 | fail |
| harmony-7 | 146 | 41 | 29 | 10 | fail |
| harmony-19 | 146 | 39 | 29 | 12 | fail |
| harmony-41 | 146 | 38 | 29 | 11 | fail |

## Every native margin failure

The following table includes all 37 supported, control-sensitive native failures.
Values are equal-donor recalls and losses; negative losses would indicate
improvement. Per-donor counts, recalls, full confusion matrices, precision/F1,
accuracy and balanced recall for every candidate remain in the archive.

| Preparation / treatment | Author label | Rare | Baseline recall | Native recall | Mean loss | Maximum donor loss |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Bcell / IFN-L1 | B intermediate | no | 0.442318 | 0.290632 | 0.151686 | 0.292466 |
| Bcell / IFN-L1 | ILC | yes | 0.371874 | 0.347985 | 0.023889 | 0.150000 |
| Bcell / IFNa | B intermediate | no | 0.542359 | 0.392921 | 0.149438 | 0.373607 |
| Bcell / IFNa | cDC2 | yes | 0.513754 | 0.513666 | 0.000088 | 0.208333 |
| Bcell / IFNb | B intermediate | no | 0.659399 | 0.586228 | 0.073171 | 0.352892 |
| Bcell / IFNb | B memory | no | 0.680633 | 0.696725 | -0.016093 | 0.157638 |
| Bcell / IFNb | HSPC | no | 0.779474 | 0.739541 | 0.039933 | 0.216535 |
| Bcell / IFNg | B intermediate | no | 0.310389 | 0.225653 | 0.084736 | 0.227564 |
| Bcell / none | B intermediate | no | 0.369336 | 0.281463 | 0.087873 | 0.213071 |
| Bcell / none | ILC | yes | 0.160193 | 0.154344 | 0.005849 | 0.165138 |
| Monocyte / IFN-L1 | CD14 Mono | no | 0.656134 | 0.723947 | -0.067814 | 0.129704 |
| Monocyte / IFNa | CD14 Mono | no | 0.712143 | 0.680337 | 0.031806 | 0.333084 |
| Monocyte / IFNb | CD14 Mono | no | 0.630129 | 0.570407 | 0.059723 | 0.357216 |
| Monocyte / IFNb | HSPC | yes | 0.734038 | 0.702198 | 0.031840 | 0.125984 |
| Monocyte / IFNg | HSPC | yes | 0.529134 | 0.516472 | 0.012662 | 0.133858 |
| Monocyte / none | CD14 Mono | no | 0.750853 | 0.755143 | -0.004290 | 0.152136 |
| NK / IFNa | CD4 TCM | no | 0.406081 | 0.365065 | 0.041016 | 0.175589 |
| NK / IFNa | CD8 TEM | yes | 0.445310 | 0.326121 | 0.119190 | 0.214286 |
| PBMC / Fresh | B intermediate | no | 0.364995 | 0.275595 | 0.089401 | 0.195890 |
| PBMC / Fresh | B memory | yes | 0.801967 | 0.777651 | 0.024316 | 0.186813 |
| PBMC / culture_IFNa | B intermediate | no | 0.845649 | 0.754719 | 0.090930 | 0.203230 |
| PBMC / culture_IFNa | B memory | no | 0.148487 | 0.290226 | -0.141740 | 0.117517 |
| PBMC / culture_IFNa | CD14 Mono | no | 0.588200 | 0.681855 | -0.093655 | 0.157844 |
| PBMC / culture_IFNa | CD16 Mono | no | 0.716291 | 0.593397 | 0.122894 | 0.359033 |
| PBMC / culture_IFNa | CD8 Naive | yes | 0.824478 | 0.778943 | 0.045535 | 0.162162 |
| PBMC / culture_no_stim | B intermediate | no | 0.429537 | 0.238104 | 0.191433 | 0.383220 |
| PBMC / culture_no_stim | CD14 Mono | no | 0.783012 | 0.687085 | 0.095927 | 0.615337 |
| PBMC / culture_no_stim | CD4 Naive | no | 0.860526 | 0.842802 | 0.017724 | 0.126545 |
| Tcell / IFN-L1 | CD4 CTL | yes | 0.821911 | 0.763057 | 0.058855 | 0.182927 |
| Tcell / IFN-L1 | CD4 Naive | no | 0.908597 | 0.867639 | 0.040958 | 0.186262 |
| Tcell / IFN-L1 | CD8 TCM | no | 0.330995 | 0.279020 | 0.051975 | 0.161716 |
| Tcell / IFNa | CD8 TEM | no | 0.440734 | 0.520811 | -0.080077 | 0.141309 |
| Tcell / IFNg | CD4 Naive | no | 0.886098 | 0.861898 | 0.024201 | 0.151610 |
| Tcell / IFNg | CD8 TCM | no | 0.372345 | 0.388809 | -0.016464 | 0.145907 |
| Tcell / none | CD4 Naive | no | 0.903016 | 0.878521 | 0.024495 | 0.108791 |
| Tcell / none | CD8 TCM | no | 0.255984 | 0.244490 | 0.011494 | 0.127329 |
| Tcell / none | NK | yes | 0.642447 | 0.640147 | 0.002300 | 0.118182 |

## Numerical verification and evidence

Five focused checks cover imbalanced class weighting, explicit held-donor
mutation invariance, direct SVD agreement, exact erasure/unseen-class behavior and
empty training. All pass. Independent augmented SVD least squares then fits
**every actual donor fold** directly from weighted per-cell rows for baseline,
native and all three Harmony matrices. Every query confusion matrix agrees
exactly. Maximum coefficient discrepancy is below 1.0e-15 and maximum intercept
discrepancy below 3.1e-14. All explicit matrix row/column coordinates and payload
hashes are checked. This verifies the diagnostic arithmetic, not label accuracy.

The original implementation computed training moments by subtracting held-donor
moments from pooled moments. The retained first attempt already passed all
per-cell checks. Before publication, training-only selection was made explicit
before summation so even floating cancellation cannot depend on query statistics.
All seven original/corrected confusion and metric results are exactly unchanged;
both executions remain archived. An earlier preparer failed before metadata
freezing because its source-hash literal was truncated; the full pinned hash was
restored and the original source remained unchanged.

This analysis is Python/NumPy and transductive: full-cohort PCA/integration precede
donor-held-out decoding. There is no new native fit, GPU benchmark, independent
study or prospective identity qualification. The native matrix is the previously
qualified seed-7 result with executable SHA-256
`e360645362ac0707cb497eaa64277e5bf792b28d79bddea707421fdcd043495d`.
All source cells, original H5AD and complete score matrices remain preserved.

The [132-member archive](evidence/2026-09-11-annotation-retention/manifest.json)
binds the exact metadata ledger, row codes, protocol, scripts, tests, all results,
all five per-cell oracles and both attempts. Manifest SHA-256:
`332d62d6f9aab85816c6e1fff8ef553a638daaa593d07c892eedf1b9dec87b08`.
Large matrices and the original source are external under explicit hashes; the
archive does not imply that those payloads are embedded. Stored and decoded bytes
are verified on both hosts.

## Reproduce

Use the existing qualified HIRISA study layout, with NumPy, h5py and ijson. Copy
the protocol to a fresh `STUDY/annotation-retention/protocol.md`, then run:

```sh
python prepare_annotation_retention.py --study STUDY --repo REPO --out STUDY/annotation-retention
python run_annotation_retention.py --study STUDY --root STUDY/annotation-retention --python PYTHON
```

The driver runs tests, seven diagnostic evaluations and five direct per-cell SVD
checks. Outputs must be new; inspect and preserve failed phases before a targeted
continuation. `archive_annotation_retention.py` packages this exact qualification,
including its retained first attempt and unchanged-result check. The general
archive verifier is `python verify_archive.py ARCHIVE`.
