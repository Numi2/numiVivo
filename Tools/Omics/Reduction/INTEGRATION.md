# Donor/batch integration and response preservation

The ridge method below remains the default. An opt-in
[native scale-aware MNN method](MNN_INTEGRATION.md) is now available through the
PCA-bundle workflow, with separate anchor witnesses and full-cohort evaluation.

The native analysis path supports transductive correction of PCA coordinates for
one explicitly selected categorical donor or batch covariate. Counts, normalized
RNA and original PCA remain unchanged. To use corrected coordinates downstream:

```json
{
  "schemaVersion": 1,
  "id": "donor-integration",
  "reduction": {},
  "integration": {"covariate": "donor", "clusters": 88, "seed": 7},
  "neighbors": {"representation": "integrated"},
  "clustering": {},
  "embedding": {}
}
```

The cluster count is explicit; 88 is the Kang benchmark setting, not an automatic
choice for every dataset. Integration requires PCA. Integrated neighbors require
integration. An omitted neighbor representation continues to use original PCA;
requesting integration alone does not silently change the graph. UMAP initializes
from the same representation used for its graph. Clustering consumes that graph.

## Algorithm and numerical evidence

The method alternates diversity-weighted soft clustering in cosine-normalized PCA
space and a mixture of weighted categorical ridge corrections. Its objective and
intercept-derived centers follow the Harmony 2 formulation inspected in the exact
`harmonypy==2.0.0` source distribution. The implementation independently uses
Double arithmetic, SplitMix64, kmeans++ and ten Lloyd initialization rounds.
It does not promise coordinate identity with the reference's Float32 and different
initialization. There are four assignment sweeps per outer iteration, in shuffled
5% blocks; probabilities use a stable log-softmax. The diversity factor for a
cluster/batch is `((2E+1)/(O+E+1))^theta`.

Each cluster fits an unpenalized intercept plus ridge-penalized batch indicators
to **original PCA**, and subtracts only the batch effects, weighted by memberships.
The intercept becomes the next cosine-normalized center. Batch levels with
cluster membership mass divided by total batch size at most 1e-5 are excluded
from that cluster's correction. Fewer than two active levels means no correction
for that cluster. The native Schur-complement solve checks normal-equation
residuals. Its complement uses a sum of positive terms to avoid cancellation.

The report retains final memberships, assignment coordinates/centers, objective
history, signed improvements, stopping reason, and ridge residual. The independent
checker reconstructs the final objective and full correction using NumPy dense
`(batches+1)`-square normal-equation solves. An objective increase is reported
separately from satisfying the relative tolerance; reaching the iteration cap is
also explicit. None of these stopping states establishes global convergence.

Missing, single-level and `unreported` covariates are rejected. A disconnected
condition/covariate design is rejected because condition and the correction
covariate cannot be separated across components. Condition labels participate in
this eligibility check and subsequent evaluation, not in fitting the correction.
Unknown options and resource-budget violations are controlled errors.

That condition-connectivity check is not a general biological-protection check.
The [full HIRISA design audit](IntegrationDesign/README.md) finds an exact overlap
between monocyte preparation and five batch labels. Donor/batch columns are also
redundant in three preparation subsets. Multiple-covariate correction remains
unimplemented; it needs explicit protected-covariate handling and preservation
qualification, not an assumption that batch effects are biologically neutral.

The preflight work index is `cells * clusters * dimensions * (maximumIterations +
10)`, capped by `maximumWork` (default 200 million). It bounds problem dimensions;
it is not a measured instruction count or performance guarantee. Dense storage
is limited to scores, memberships, centers and small batch systems, never a
cells-by-genes expression matrix. The original analysis route remains resident.
[File-backed integration](FILE_INTEGRATION.md) now shares the same solver with
bounded latent-matrix mappings and explicit integrated graph inputs. Both use
one categorical correction covariate, at most 100 clusters and 128 levels;
neither is an unseen-donor transform. An explicit work maximum up to 100 billion
is supported, with the unchanged 200-million default.

## Experimental benchmark

`check_integration.py` uses all 2,651 retained Kang B cells from eight paired
donors, with control and interferon-stimulated samples. It requires every donor
to contain both conditions and verifies cell identities against the retained RNA
matrix. Harmony sees only donor labels and the same native 20-PC representation.
Reference seeds are 7, 19 and 41, one CPU thread, fixed theta 2, ridge 1,
temperature .1, 88 clusters and at most ten outer iterations. Exact parameters,
package versions, source report hash, coordinate hashes and objective histories
are retained. The reference dependency is external and optional for native use.

Evaluation uses exact 30-neighbor searches in 64-row distance blocks:

- Donor mixing is computed **within each condition**, using normalized donor
  entropy and excess same-donor neighbor fraction over the actual eligible donor
  composition. Metrics average donors equally. Negative excess is reported as
  such; it is not evidence of optimal or biologically correct mixing.
- A standardized logistic classifier trains on seven donors and tests the eighth.
  Scaling is fitted on training cells. Each fold and balanced accuracy are retained.
  Integration is transductive and saw the unlabeled test cells, so this is not a
  held-out perturbation or prospective unseen-donor prediction benchmark.
- An explicitly declared interferon program averages log-normalized IFI6, IFIT1,
  ISG15, MX1 and ISG20. No missing gene is silently dropped. Neighbor-smoothed
  program correlation is measured per donor, and also within donor/condition
  strata so condition separation cannot be its only support.
- A deliberately label-informed condition-centering control removes the mean
  response in PCA. Its unchanged within-condition distances expose the inability
  of mixing metrics alone to detect response erasure. It must reduce classifier
  accuracy by at least .10. It is an evaluation control, never an integration option.

The declared cohort gates require reduced same-donor excess, no more than .02
loss of balanced condition accuracy, and no more than .05 loss of either program
correlation. Raw values and failed gates are retained. These margins are engineering
checks for this cohort, not statistical confidence intervals or general biological
acceptance. PBMC3k has one library and is ineligible for donor integration. Diverse
cell types, rare states and additional independent donor-resolved studies still
need qualification before any competitive multi-donor claim.

```sh
python check_cli.py --integration --clustering --embedding --binary /path/to/numivivo --imported /path/to/kang/imported --out /new/native
python check_integration.py --report /new/native/report.json --out /new/reference.json
python check_neighbors.py --report /new/native/report.json --out /new/neighbors.json
```

Primary references: [Harmony paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6884693/),
[harmonypy](https://github.com/slowkow/harmonypy), and the exact 2.0.0 source
archive recorded in the qualification evidence. See `evidence/2026-09-09-integration/`.

## Qualified production run, 2026-09-09

All 35 native tests in nine suites passed, as did 11 production workflow checks
(including replay, exact repeated receipt and controlled failure cases) and the
existing 18 CLI assertions across 24 commands. Original processed data and PCA
are exactly equal to the previous unintegrated production report. Binary identity
is `d696896b1c87b0eab2c71923928015234e2f0cb049d259b08ca6c98f876ff6b8`;
source hashes, the initial fixture compilation failure and final logs are retained.

| Donor-balanced measurement | Original PCA | Native, seed 7 | Harmony reference range, three seeds |
|---|---:|---:|---:|
| Within-condition donor entropy | 0.801481 | 0.860630 | 0.858571–0.861234 |
| Excess same-donor fraction | 0.027334 | -0.015533 | -0.015490–-0.014152 |
| Cross-donor condition balanced accuracy | 0.991984 | 0.994469 | 0.991291–0.993777 |
| RNA program neighbor Spearman | 0.895781 | 0.888251 | 0.887094–0.889888 |
| Within-condition program Spearman | 0.584323 | 0.551650 | 0.542384–0.556746 |

Both program correlations decrease; they pass the declared tolerance, not a claim
of complete preservation. The response-erasure control leaves within-condition
mixing unchanged but lowers cross-donor classification accuracy to 0.353368.
The native objective reached its relative tolerance after three corrections;
independent objective reconstruction gives 632.191497888315 versus native
632.191497888196. Independent ridge correction has maximum coordinate error
2.176e-14. This older B-cell run tested native seed 7 only. The newer
[full-cohort evaluation](FILE_INTEGRATION.md) tests three native and reference
seeds on full Kang and independent Hagai, retaining the Kang NK-cell failure.

The integrated 15-slot neighbor graph has exact independently verified membership
and ordering, 55,508 fuzzy connectivity entries and one connected component.
Its maximum reference weight error is 4.457e-6. Downstream Louvain has nine
communities, modularity 0.640633, and no disconnected communities; this is not
cell annotation. UMAP trustworthiness is 0.917899 against 0.917036–0.918928 for
three reference layouts; neighbor recall is 0.233094. Most original neighbors
are not retained in two dimensions. These graph/layout checks use integrated
coordinates explicitly and do not add biological qualification.
