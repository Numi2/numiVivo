# Arc Virtual Cell Challenge 2026

Native NumiVivo input/submission integration plus a separate, pinned official
scoring adapter. This is a benchmark implementation, **not a trained zero-shot
model, a prediction result, or an uploaded competition submission**.

## Frozen specification

The adapter targets `ArcInstitute/cell-eval2` commit
`5e64833518a6603a0301cbe28185d49c30f4a986`, package **0.16.0**, competition **rule 3**,
and its complete `vcc2026` preset. See `evaluator-lock.json`.

Official sources:

- https://arcinstitute.org/news/virtual-cell-challenge-2026
- https://github.com/ArcInstitute/cell-eval2/blob/5e64833518a6603a0301cbe28185d49c30f4a986/docs/vcc2026_metrics/vcc2026-metrics-brief.md
- https://github.com/ArcInstitute/cell-eval2/blob/5e64833518a6603a0301cbe28185d49c30f4a986/src/cell_eval2/configs/vcc2026.yaml
- https://github.com/ArcInstitute/cell-eval2/blob/5e64833518a6603a0301cbe28185d49c30f4a986/CHANGELOG.md

Each phase has three contexts, each with 300 targets and the ordered common
18,533-gene axis. Predict from basal non-targeting controls and target identities,
not perturbed outcomes from those contexts. Submissions contain raw nonnegative
integer counts in `X`, all targets plus `non-targeting`, and no row above
1,000,000 counts. **Predicted cell counts need not equal 400**, and predicted
library depth is not required to equal 20,000: those are reference preparation
quantities. Native admission never normalizes, rounds, fills genes, intersects
axes or repairs missing target predictions.

Arc's six scored metrics are `pds_cosine`, `expr_mse_unbiased_capped_norm`,
`de_wilcoxon_direction_fidelity_yield_raw`, `de_wilcoxon_direction_reach_raw`,
`de_wilcoxon_sig_jaccard`, and `de_wilcoxon_lfc_nmae`. The preset also emits four
unscored diagnostics. PDS excludes the entire target panel; other members exclude
their own target. Scoring uses reference controls, not predicted controls, and
the preset owns the Wilcoxon/filter/BH, 50,000 bulk normalization, sampling
correction, ties, replicate anchors and per-metric clipping rules.

No copies of those metric algorithms are introduced here. The official scorer's
`from_replicate` column is authoritative; `from_baseline` at `avg_score` must be
null. All six metric scores and all three context scores are required for a phase
mean. Incomplete or failed execution cannot emit a complete phase score.

## Ownership

`VivoArc2026Contract` owns the versioned native contracts and raw-count bounds.
`VivoArc2026H5AD` reuses `VivoOmicsFileSnapshot`, `VivoH5ADFrameReader`,
`VivoH5ADCountReader` and `VivoCanonicalJSON`. Existing CSR/CSC/dense count decoding
is retained; the matrix is not materialized. Cell IDs, labels and row totals are
resident. Current local limits are 2,000,000 cells, 2,000,000,000 nonzeros and
64 GiB/file, not additional Arc competition limits. Unsupported AnnData encodings
fail explicitly. Native exact-integer validation is intentionally stricter than
upstream floating-point integer-tolerance admission.

`arc2026-*` commands are thin public CLI routes. `evaluate.py` is offline reference
orchestration and never participates in native inference or biological simulation.
No second model, count runtime, training pipeline or dependency on Python in the
production neural/physical loop is added.

## Native preparation and submission

Use the actual released context identifiers and files. A/B/C in the templates
are placeholders, not claims about official cell-line identities. Paths resolve
relative to the invocation's working directory; absolute paths are recommended.
`targetsJSON` is the adapter's explicit array of 300 exact target label strings.
When translating an official target list, preserve the original release and
record that conversion in source provenance; do not infer targets from RNA.

```sh
numivivo arc2026-prepare prepare.json --output /new/query
numivivo arc2026-verify-query /new/query
# Run your predictor separately using ONLY /new/query and reviewed training inputs.
numivivo arc2026-pack /new/query --plan pack.json --output /new/submission
numivivo arc2026-verify /new/submission
```

A query contains unchanged control H5ADs, exact target lists, hashes, common gene
identities and independently reconstructable count summaries. It accepts no
reference-outcome or scoring-bundle path. A submission includes the frozen query,
three `contexts/<id>/prediction.h5ad` files and a source-bound model/training
provenance declaration. It publishes only after every context validates. Existing
outputs are never overwritten; verification re-reads private immutable snapshots.

The model artifact hash and excluded-context list are **declarations**, not proof
that particular weights generated the predictions or that all historical training
was uncontaminated. The training manifest is retained verbatim. Review its data
sources, use of validation feedback and context exclusions. Native checks do not
forensically audit arbitrary AnnData auxiliary metadata. These receipts establish
supported-file input conformance, not competition eligibility or biological skill.

## Exact offline scoring

Use a Python environment with the pinned official evaluator installed. The base
package supports CPU scoring; additional engines must match the supplied bundle.
For example, prepare a separate checkout/environment without changing NumiVivo:

```sh
git clone https://github.com/ArcInstitute/cell-eval2.git /path/cell-eval2
git -C /path/cell-eval2 checkout 5e64833518a6603a0301cbe28185d49c30f4a986
python3 -m venv /path/arc-eval-env
/path/arc-eval-env/bin/pip install /path/cell-eval2
# Install the bundle's DE backend/extras as specified by Arc, when required.
/path/arc-eval-env/bin/python Tools/Omics/Arc2026/evaluate.py run \
  --native /path/to/numivivo --evaluator /path/cell-eval2 \
  --submission /new/submission --references references.json --output /new/evaluation
```

The reference plan is separate from prediction. It explicitly binds each query's
control/target hashes to a reference-file hash and a bundle-manifest hash. Fill
those bindings from reviewed release provenance before running; filenames and
matching dimensions are not evidence of context identity. All bundle files are
then snapshotted and hashed, not just the upstream manifest. This does not turn a
user-supplied bundle into an authenticated official download.

The adapter copies and verifies all submitted predictions before reading outcomes.
It invokes native submission reconstruction, checks the evaluator's actual Git
blob contents against the pinned commit, and runs only a frozen copy of those
tracked sources (not checkout bytecode). It records dependency versions, native
binary and adapter identities, complete logs, reference/prediction/bundle hashes,
raw metric outputs and official scaled scores.

It runs `cell-eval2 run --preset vcc2026`, then `cell-eval2 score --real-bundle`,
without supplied DE tables, metric overrides, mismatch waivers, altered anchors or
fallback scales. Device and DE backend are selected from the bundle; Arc's exact
submission/bundle compatibility checks remain mandatory. A 0.16.0 evaluator must
reject old 0.15.0 bundles even though some prose describes those historical
bundles. **Never restamp a bundle or silently replace its GPU engine with a CPU
engine.** Backend availability or compatibility failures retain logs and no score.

If reference outcomes are withheld, prepare the native submission and stop before
local scoring. This adapter does not retrieve hidden outcomes, invent reference
bundles, register a team, automate portal submission, or implement an undocumented
portal archive convention. It produces validatable per-context H5AD files; the
competition portal remains responsible for official eligibility and acceptance.
A local score is not a leaderboard score. Fresh references and bundles are still
needed to qualify real evaluation here.

The reference plan may list fewer than three available scoring contexts. Individual
results can be retained, but the phase score remains null. The adapter never averages
validation and final-test contexts together. Runs use create-only output directories;
failed attempts and logs remain available. Re-run into a new directory rather than
mutating an existing evaluation.

```sh
python Tools/Omics/Arc2026/evaluate.py verify /new/evaluation
```

This verifies artifact hashes and recomputes six-metric/three-context aggregation;
it does not rerun Wilcoxon, the metric kernels or native HDF5 inspection. A fresh
`run` is required for numerical evaluator reproduction.

## Checks and current execution boundary

```sh
python -m unittest discover -s Tools/Omics/Arc2026 -p test_adapter.py -v
# Supported Apple host, after building the current numivivo executable:
swift test --filter VivoArc2026Tests
python Tools/Omics/Arc2026/check_native.py --native /path/to/numivivo --output /new/arc-native-check
```

The native executable check creates explicitly synthetic 18,533-gene/300-target
H5ADs. It tests complete phase preparation/pack/replay, non-400 cell counts,
fractional counts, excessive depth, axis changes, missing targets, immutable output
and altered prediction bytes. It does not use real outcomes or produce a score.

Development in this increment executed 15 standard-library adapter tests and Swift
syntax parsing. Full Apple type checking/HDF5 execution, the authored native tests,
installation/execution of the official evaluator and scoring on Arc's released
contexts were **not run**. There is no retained benchmark score or performance claim.
