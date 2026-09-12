# Variance-retaining RNA profile: preservation gates fail

The native candidate retains RNA PCA scores (UΣ), normalizes each cell row to unit length, and leaves ATAC preprocessing and all 30 components per modality unchanged. It uses the qualified native neighbor and modality-weight kernels on all 2,711 nuclei. The production default remains unchanged.

| Candidate weighted reconstruction versus baseline | Neighbor-mean SSE change | Fuzzy-mean SSE change |
| --- | ---: | ---: |
| Original RNA-only | +0.64676% | -0.36732% |
| Matched variance-retaining RNA-only | +0.70410% | +0.06928% |
| Original/matched ATAC-only | +0.08414% | -0.33287% |

Positive means worse. The frozen rule requires no worsening in aggregate SSE, macro variance-normalized error and every nonempty detection-frequency stratum, for both operators. RNA must meet both original and matched baselines; ATAC must meet the original baseline. Both modality gates fail. Earlier [RNA](../RNAPreservation/README.md) and [ATAC](../ATACPreservation/README.md) failures remain valid.

RNA evaluates all 34,601 genes outside the selected 2,000 HVGs. ATAC evaluates all 98,319 peaks. The candidate RNA-only baseline itself improves aggregate reconstruction, so comparison only with the original RNA baseline would hide a failure against the matched baseline. Equal-weight joint results, all per-feature errors, null errors and support groups are retained.

## Evidence and limits

The native harness completed successfully. Independent SciPy exact neighbor IDs and distances pass for RNA and ATAC. Independent sparse shared-neighbor formulas agree with native modality weights within 2.03e-15 (score error 1.85e-13). Candidate joint/weighted-neighbor and fuzzy-reference comparisons have not been rerun; this is partial numerical verification, not complete candidate qualification.

Squared-error calculations agree with a second algebraic calculation and scalar checks on 17 prespecified features per modality. RNA maximum algebra/scalar differences are 5.37e-10/4.55e-12; ATAC differences are 7.17e-12/1.94e-12.

The manifest binds the archive containing the pre-execution protocol, native harness and graphs, reference code/report, scoring scripts/protocols/reports, per-feature arrays and combined decision. Original inputs and qualified library hashes are bound by the protocols. Archive SHA-256 is checked on publication.

These are development reconstruction diagnostics on one already inspected donor. HVG selection and normalization used all RNA genes; LSI used all ATAC peaks. They do not establish independent held-out prediction, regulatory causality or clinical outcomes. General biological prediction remains unproven. This candidate does not justify changing the product default or relaxing either gate. Further progress requires a separately specified experimental evaluation rather than choosing a favorable graph/operator after inspecting these outcomes.
