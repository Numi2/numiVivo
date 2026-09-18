# Published binder baseline and native structural increment — 2026-09-18

## Measured conclusion

**The fixed three-score ensemble does not improve top-10 binder selection over
the fixed Boltz-2 ipSAE baseline in this selected retrospective panel.** It ties
four target/assay comparisons and loses two. The ensemble is not promoted.

This is an actual published-table result, not a synthetic benchmark. It establishes
a comparison that future Numi structural/physical analysis must improve; no physical
features enter these six published evaluations and no physical uplift is demonstrated.

## Data and frozen protocol

Original source: Anthropic `claude-protein-binder-design`,
revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`,
`data/tables/design_summary.csv`, CC BY 4.0.

- Original table: 1,440 rows, 1,213,035 bytes; the complete source was retained.
- Explicitly selected panel: 270 designs across BBF-14, MBP and EGFR.
- Binary evaluable observations: 256 for Adaptyv and 267 for Twist; these are
  separate measurements of overlapping candidates, not 523 independent designs.
- Source Git blob: `d1573ba03e8322c70ccb3a40e46418e86b40e2dc`.
- Source SHA-256: `a96dc3a4fcc293dd393a8c56fa367dadcac4db326c30f1c3ea3d8c2b9f72367a`.

The six fixed plans hold out one complete target at a time for each assay. The
baseline is `ipsae_min_boltz2`; the ridge-logistic ensemble adds
`ipsae_min_ptxv2` and `ipsae_min_ef2full`. The ridge penalty is 0.1 and K is 10.
Normalization and fitting use training rows only. Exact-sequence groups shared
with any held-out row are purged from training; unavailable outcomes are not
converted into negatives. There was no parameter search in this campaign.

The runner writes all plans before import/scoring. The [protocol](protocol.json)
and plan hashes bind the executed configuration; they do not prove that no human
has previously inspected public outcomes. These targets are now reused development
data for subsequent feature experiments, not a fresh external validation set.

## All six results

Top-10 hits mean experimentally classified binders among ten selected candidates,
with fractional selection at score ties. The values happen to be integers here.

| Assay | Held-out target | Training N | Test N | Test binders | Boltz-2 hits@10 | Ensemble hits@10 | Boltz-2 AUROC | Ensemble AUROC |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Adaptyv | BBF-14 | 178 | 78 | 1 | 0 | 0 | 0.090909 | 0.155844 |
| Adaptyv | MBP | 168 | 88 | 0 | 0 | 0 | unavailable | unavailable |
| Adaptyv | EGFR | 166 | 90 | 10 | 2 | 0 | 0.585000 | 0.320000 |
| Twist | BBF-14 | 178 | 89 | 3 | 1 | 0 | 0.399225 | 0.441860 |
| Twist | MBP | 179 | 88 | 1 | 0 | 0 | 0.655172 | 0.747126 |
| Twist | EGFR | 177 | 90 | 10 | 2 | 2 | 0.595000 | 0.640000 |

The [complete summary](summary.json) also retains average precision, Brier scores
and the training-prevalence control. In particular, a higher AUROC in some folds
does not erase the lack of improved top-10 selection. The Adaptyv MBP endpoint has
no measured binders, so it cannot evaluate discrimination. Several other folds
have only one or three positives. No confidence interval, clinical claim or broad
generalization follows from these counts, and the six comparisons are not treated
as six independent biological replications.

## Actual execution and replay

Initial campaign commit: `8c125e4a83b61be2b37a6b7bbc7f773b54412823`.
Initial successful CI run: https://github.com/Numi2/numiVivo/actions/runs/35399912744

Current tested implementation: `cbcf613b5a283f2a45b28ef010006195f5b0e984`.
Successful focused CI run: https://github.com/Numi2/numiVivo/actions/runs/35401015828
Job: `105780669224`. Swift 6.2.1 on x86_64 Linux, Python 3.12.3.

The current CI completed 50 Swift tests with zero failures, seven Python tests,
the original 90-row synthetic score campaign, a 17-structure synthetic campaign
with matched 12-training/five-test rows, and six mocked-source campaign folds.
It then downloaded the actual pinned table and executed all six real folds.
Python independently checked the imported source records, eligible populations,
scaling, model gradients, predictions and ranking/probability metrics. Native
bundle verification reconstructed every import and evaluation.

Current CI artifact: https://github.com/Numi2/numiVivo/actions/runs/35401015828/artifacts/10571040764
ZIP SHA-256: `fb62f26cb415f03f81a39a2e5157dd44bed7166e681dc434a5ac8d1c9ceb9a47`.
The downloaded ZIP and all 85 entries in its manifest were verified. The artifact
contains 86 files including its manifest. CI artifact retention is 14 days;
this committed summary/protocol is permanent source history, not a substitute
for preserving the full campaign archive when long-term replay is required.

A separately built local executable reran the current campaign from the verified
source bytes. Every complete JSON evaluation report matched the current CI report
exactly, including training models and held-out predictions, not just rounded
summary metrics. The local and CI executables have distinct identities, recorded
in [execution.json](execution.json). Current and initial CI fold metrics also agree.

## Scientific boundary and next increment

The new structural code computes static geometry from the existing molecular
container and compares a geometry-augmented model with a score-only control on
identical candidates. It passed synthetic independent geometry/fit checks, but
actual published structures have not yet passed through that route. Target/model
origin remains caller-declared, and the resident archive has explicit work/memory
bounds. See [STRUCTURAL_ANALYSIS.md](../../STRUCTURAL_ANALYSIS.md).

No Boltz-2 weights were run, no MD-derived features were generated, no complete
Apple/Metal application was built, and no measured improvement in candidate selection
was established. Focused binder CI success is not a claim about every repository
workflow. The next increment is source-bound actual structure import and a frozen
matched comparison, followed by independently qualified preparation and sampling.
