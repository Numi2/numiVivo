# Outcome-blind fixed ranking

`VivoBinderRanking` separates fixed ranking from experimental assessment. The existing
`binder-evaluate` fitted-model reports and outcome-conditioned cohort semantics are unchanged.

```sh
numivivo binder-ranking-query IMPORT_BUNDLE NEW_QUERY
numivivo binder-rank NEW_QUERY Tools/BinderBenchmark/plans/ranking-confidence3.json NEW_RANKING
numivivo binder-assess-ranking IMPORT_BUNDLE NEW_RANKING NEW_ASSESSMENT
numivivo binder-verify NEW_ASSESSMENT
```

All bundles use the existing `VivoBinderBundleIO` publisher, canonical JSON, SHA-256,
no-overwrite and same-executable replay. Query/ranking bundles contain no raw CSV,
assay, experimental outcomes or preserved source fields. Unknown fields in query/plan
JSON are rejected. This is typed separation, not an operating-system sandbox or proof
that a human never inspected the public outcomes. A standalone query's claimed source
digest is authenticated against the original CSV only when assessment reconstructs it.

`fixedScore` retains a named upstream normalized score. `targetStandardizedMean`
computes a positive-weighted mean of population-SD z scores separately within each
target's complete-case query pool. Constant components contribute zero; a target with
no complete candidates is recorded explicitly. All outcomes, including unavailable
ones, remain in the query pool. Normalization is transductive/pool-dependent, not a
training-derived calibration or binding probability. Missing predictors are not zeros.

The committed confidence3 plan uses equal weights for `ipsae_min_ef2fast`,
`ipsae_min_ef2full`, and `ipsae_min_ptxv2`. This is the confidence-only profile motivated
by Anthropic's report, applied to the published seed-best re-scoring table. It does not
reproduce campaign-time inputs/inference or optional self-consistency DockQ terms, and
does not measure equal-compute inference performance. Arbitrary weights cannot be
searched on held-out outcomes and then represented as a fixed reference.

Top-K is chosen before outcome joining. Exact score ties receive fractional inclusion
weights; no identifier-based tie winner is presented as evidence. Assessment preserves
selection weight by every outcome category, measured binder hits, nonbinary weight,
and precision conditional on selected binary measurements. It never removes an
untested selected candidate and replaces it with a lower-ranked tested one. Binary
AUROC/AP use all ranked candidates with binary outcomes, without changing selection.
These semantics differ deliberately from the older outcome-conditioned benchmark;
compare methods on a common pool, not across those different endpoint definitions.

## Checks and actual-data runner

```sh
python3 Tools/BinderBenchmark/native_package.py /tmp/binder-ranking-native --test
python3 Tools/BinderBenchmark/check_rankings.py /tmp/binder-ranking-native/.build/debug/BinderCLI
python3 Tools/BinderBenchmark/run_rankings.py SOURCE.csv \
  /tmp/binder-ranking-native/.build/debug/BinderCLI /tmp/binder-ranking-campaign
```

The actual-data runner requires the existing pinned source SHA/Git-blob identity,
retains all source/plan/result files, evaluates both fixed plans against the separate
Adaptyv and Twist endpoints, checks native normalization/ranking/metrics independently
in Python, and requires outcome-blind query/ranking bytes to agree across laboratories.
It is restricted to the existing 270-design BBF-14/MBP/EGFR panel. The original
1,440-row CSV is preserved; no other target is modelled or selected by this runner.

This work does not add structural/MD features, fit-support calibration, a trained
rank correction, model-weight inference, or an experimentally validated improvement.
Portable tests compile exact binder source but do not qualify the full Apple/Metal
product. Evidence remains tied to the checked source revision and executable.

Source documentation: https://huggingface.co/datasets/Anthropic/claude-protein-binder-design
and `structure_and_pae/README.md` at revision
`9e1b81696da46835e9e9cde9a3da976e0abc92ab`. Dataset attribution: Anthropic, CC BY 4.0.
