# Native sparse UMAP-compatible coordinates

Add an explicit embedding request alongside reduction and neighbors:

```json
"embedding": {
  "dimensions": 2,
  "epochs": 500,
  "minimumDistance": 0.1,
  "spread": 1,
  "learningRate": 1,
  "negativeSampleRate": 5,
  "repulsionStrength": 1,
  "seed": 7,
  "maximumUpdates": 200000000
}
```

The existing native analysis, export and reconstruction route carries the result.
Embedding does not replace the PCA representation or the graph used for clustering.
It requires explicit PCA and neighbor options. Older plans omit it by default.

The differentiable UMAP distance curve `1 / (1 + a*d^(2*b))` is fitted natively to
300 samples of the offset exponential target. Log-positive parameters and a damped
Gauss–Newton solve require a stationary fit; there is no fallback to a guessed
curve. The first requested PCA dimensions are scaled independently to [0,10] and
jittered uniformly by at most 1e-4, using the declared seed. These initial
coordinates are retained alongside the final coordinates and curve parameters.
Two or three dimensions are supported.

Optimization follows the Euclidean UMAP attraction/negative-sampling equations,
including coefficient clipping at ±4, the 0.001 repulsive-distance regularizer,
updates to both endpoints of attractive edges, and head-only negative updates.
Edge sample periods derive from their weights. Edges below max-weight/epochs
are omitted from optimization and their count is reported. The learning-rate
schedule and 20–2000 fixed epoch range are explicit. This implementation uses
Double coordinates and SplitMix64 sampling, so it does not promise bit-identical
trajectories to umap-learn's Float32/Tausworthe implementation.

The graph and schedules remain sparse; only cells-by-2/3 coordinate arrays are
dense. The configured work cap bounds both edge visits and the estimated/actual
attractive-plus-negative update counts. Exhausting it rejects publication rather
than silently reducing the epochs. Coordinates must remain finite. The result
records update counts, retained/discarded edges, initial and final coordinates,
cell identities and maximum coordinate magnitude. Fixed epochs are not a claim
of optimizer convergence or a monotonic full cross-entropy objective.

## Evidence boundaries

`check_embedding.py` compares the curve against SciPy's UMAP reference fit and
independently reconstructs the full sample schedule. On the exact same graph and
initial coordinates it runs three umap-learn optimizer seeds (7,19,41). Quality is
measured at 15 neighbors against the original 20-dimensional PCA coordinates,
using trustworthiness, continuity and neighborhood recall. Both reference and
quality calculations avoid a full dense cells-by-cells matrix; distance ranks
are computed in blocks of 64 cells.

Predeclared engineering gates require native trustworthiness at least 0.90 and
within 0.02 of the lowest reference value, and recall within 0.05 of the lowest
reference value. All seed results and initialization metrics are retained. These
are comparisons on the tested datasets, not biological validation or automatic
acceptance for other data. A visually separated group is not a cell-type label;
UMAP distances, apparent densities and cluster sizes are not calibrated physical
or biological measurements.

```sh
python check_cli.py --clustering --embedding --binary /path/to/numivivo --imported /path/to/imported --out /new/native
python check_embedding.py --report /new/native/report.json --out /new/reference.json
```

The initial scope is all 2,651 Kang B cells and all 2,700 PBMC3k cells from the
previous count/graph qualifications. The analysis still inherits resident count
limits and the exact-neighbor pair budget. Spectral initialization, transformation
of unseen cells, densMAP, integration and out-of-core execution remain open.
Biological preservation and cross-seed/donor robustness need further evaluation.
Use stored cell identities when mapping coordinates back to AnnData axes.

Final production evidence is under `evidence/2026-09-09-embedding`: 32 Swift tests
in eight suites, nine workflow assertions per dataset, and 18 existing CLI
assertions across 24 commands. Both real-data reference gates passed.

| Dataset | Native trustworthiness | Reference range | Native 15-neighbor recall |
| --- | ---: | ---: | ---: |
| Kang B cells | 0.9011 | 0.9023–0.9032 | 0.1984 |
| PBMC3k | 0.9395 | 0.9397–0.9409 | 0.2256 |

The low neighbor recalls remain visible: two-dimensional coordinates do not
retain most original 15-neighbor relationships. The full reference reports also
retain continuity, initialization metrics, every seed and exact update counts.
[Rendered coordinates](evidence/2026-09-09-embedding/coordinates.png) use declared
condition for Kang and numeric graph-community IDs for PBMC3k; no cell types are
inferred. The figure is linked to source report hashes in `figure-source.json`.

`publication-scope.json` records preserved upstream adapter/documentation changes
that do not alter the tested native compilation inputs. Raw compiler logs retain
existing warnings and original whitespace.
