# Outcome-blind rankings and training-support gate — 2026-09-18

Execution occurred on September 18 in America/Los_Angeles (September 19 UTC).
This record concerns ranking/support development, not improved binder selection.

## Delivered code

- `52f387ff`: outcome-blind fixed ranking and separate experimental assessment.
- `b156046b`: public CLI, source-replayable ranking bundles, fixed plans and campaign.
- `58556d17`: opt-in training-group/target support policy using the exact legacy cohort owner.
- `17d08e03`: supported-fit CLI, immutable refusal reports and independent support checker.
- `a9584261`, `2a0e954b`: focused CI coverage for real rankings and explicit support decisions.

Concurrent raw PDB/mmCIF source transport, parsing and CLI evaluation were preserved
and tested with this work. They are not reattributed to the ranking implementation.
Non-fast-forward writes were rejected and reconciled without force pushing.

## Executed checks

The combined native source through `1c24261e8e74c07b19f79d5f8079c7b99058caa1`
passed **90 Swift tests, zero failures**, and 17 offline Python tests (seven original
source/plan tests plus ten source-transport tests). Exact checked source hashes,
toolchain and executable identity are in [execution.json](execution.json).

Additional actual-file checks passed: the 90-row synthetic fitted CLI campaign;
17 synthetic structures through embedded and original PDB/mmCIF paths; independent
Python geometry/fitting/metric checks; the outcome-blind ranking CLI campaign;
all six real legacy folds; 12 real ranking assessments; and six real support decisions.

The latest combined local binary has SHA-256
`c62488282b873ddf48a886a6b01f6b620877f91afe9f2589a6a4a33d67834c8e`.
All six complete legacy `report.json` files remained byte-identical to the original
`cbcf613b` CI campaign after the shared-cohort refactor. Accepted support-gated
fitted reports also matched their original legacy reports exactly.

These checks compile the actual binder/structure owners and CLI with the portable
canonical-JSON/OpenSSL facade and native artifact primitives. They do not build
the complete Apple/Metal product or exercise its production artifact store.

Focused CI at `a9584261` completed successfully in run `35415244187`, covering
ranking and raw-source integration before the support CLI. The later own CI run
`35415908432` at `2a0e954b` was cancelled; no completed support-CI result is inferred
from it. The final combined support results recorded here are local execution
results. Future CI completion does not retrospectively change this execution record.

## Actual outcome-blind reference rankings

Source: Anthropic `claude-protein-binder-design`, CC BY 4.0, revision
`9e1b81696da46835e9e9cde9a3da976e0abc92ab`. The original 1,440-row table was retained;
this campaign selects the existing 270 designs across BBF-14, MBP and EGFR.
Source SHA-256: `a96dc3a4fcc293dd393a8c56fa367dadcac4db326c30f1c3ea3d8c2b9f72367a`.

All 90 candidates per target enter fixed ranking before experimental outcomes
are joined. The confidence3 reference uses equal-weight within-target population-SD
standardization of published seed-best ESMFold2-fast/full and Protenix ipSAE scores.
It is confidence-only re-scoring, not new model inference, optional self-consistency
terms, exact campaign reproduction or equal-compute comparison.

Query and ranking bytes were **identical across the two laboratory endpoints**.
Missing outcomes do not change normalization or cause top-ten selections to be
refilled. Cutoff ties receive equal fractional inclusion rather than an ID winner.
The comparison is operational selection, not the older binary-outcome-conditioned
benchmark; do not mix those cohort definitions.

| Assay | Target | Fixed Boltz-2 known hits@10 | Confidence3 known hits@10 | Boltz-2 selected nonbinary weight | Confidence3 selected nonbinary weight |
| --- | --- | ---: | ---: | ---: | ---: |
| Adaptyv | BBF-14 | 0 | 0 | 4 | 0 |
| Adaptyv | EGFR | 2 | 2 | 0 | 0 |
| Adaptyv | MBP | 0 | 0 | 0 | 0 |
| Twist | BBF-14 | 1 | 0 | 0 | 0 |
| Twist | EGFR | 2 | 2 | 0 | 0 |
| Twist | MBP | 0 | 0 | 0 | 0 |

**No improved top-ten binder selection was established.** A nonbinary selected
outcome is not a demonstrated non-binder. Multiple laboratory measurements and
multiple predictor samples are not independent experimental designs. The already
inspected panel remains development data, not an independent test set.

## Explicit support-policy execution

The committed `support-coverage-example.json` is an illustrative policy requiring
at least two exclusively positive and two exclusively negative declared training
groups, two training targets with both group classes, and at most three features.
It is not a calibrated power, accuracy or minimum-sample-size rule.

| Held-out fold | Eligible training candidates | Positive groups | Example policy decision |
| --- | ---: | ---: | --- |
| Adaptyv BBF-14 | 178 | 10 | Decline: both-class target coverage missing |
| Adaptyv EGFR | 166 | 1 | Decline: positive-group and target coverage missing |
| Adaptyv MBP | 168 | 11 | Experimental fit admitted |
| Twist BBF-14 | 178 | 11 | Experimental fit admitted |
| Twist EGFR | 177 | 4 | Experimental fit admitted |
| Twist MBP | 179 | 13 | Experimental fit admitted |

Declined runs publish complete counts/reasons and exit 2 without a fitted model.
The separately retained fixed ranking is not modified. Accepted runs execute the
unchanged legacy fitter; numerical errors are not hidden as successful fallbacks.
Policy decisions and counts were independently reconstructed in Python. Rehashing
a changed policy without recomputing its dependent report is rejected by replay.

Admission is not model promotion. Coefficient stability, calibrated uncertainty,
near-homology/design-family independence and external validation remain open.
The gate is opt-in; historical `binder-evaluate` behaviour is deliberately unchanged.

## Reproduce

```sh
python3 Tools/BinderBenchmark/native_package.py /tmp/binder-native --test
BIN=/tmp/binder-native/.build/debug/BinderCLI
python3 Tools/BinderBenchmark/check_rankings.py "$BIN"
python3 Tools/BinderBenchmark/fetch_source.py /tmp/binder-source
python3 Tools/BinderBenchmark/run_published.py /tmp/binder-source "$BIN" /tmp/binder-legacy
python3 Tools/BinderBenchmark/run_rankings.py /tmp/binder-source/source.csv "$BIN" /tmp/binder-rankings
python3 Tools/BinderBenchmark/check_support.py /tmp/binder-legacy "$BIN" \
  Tools/BinderBenchmark/plans/support-coverage-example.json /tmp/binder-support
```

Use new output directories. All bundles require their producing executable for
strict replay; rebuilding creates a new identity and requires fresh source import.
The separately retained replay archive is 11,450,129 bytes with SHA-256
`52938b4435533bfc46127fa6d409d3e0216c07c9e7be28032dd341489cb2af1e`.
It includes original source, plans, individual reports, executable and check logs;
this committed summary is not a substitute for the complete replay artifacts.

## Next scientific boundary

This increment adds no predictor-weight inference, multimodel PAE/seed feature
aggregation, new burial/packing observations, learned structural correction,
MD-derived feature or independent biological result. The next selection experiment
must add source-bound information beyond these fixed scores, preserve the candidate
pool, and compare against both fixed references. Concurrent raw-structure work can
supply inputs through the same existing owners; it does not itself establish an
improved ranking. Full macOS/Metal qualification remains separate.
