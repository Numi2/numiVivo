# Binder workflow: executed portable check, 2026-09-18

## Result

**35 Swift tests passed, zero failures. Seven offline Python tool tests passed,
zero failures. The complete synthetic CLI campaign passed.**

Environment: Swift 6.2.1 on x86_64 Linux. The checker builds the exact native
`NumiVivoKit/Binder` owners and `VivoBinderCLICommands` in a temporary Swift package.
The existing `Tools/Posterior/PortableSupport.swift` supplies the portable
OpenSSL SHA-256/canonical-JSON facade. The full Apple product and its production
artifact-store implementation were not built or qualified.

The 90-row synthetic campaign uses 60 training candidates and 30 candidates from
a wholly held-out synthetic target, three score features and seed 183. It runs
`binder-import`, `binder-evaluate` and `binder-verify` against actual files.
An independent standard-library Python calculation checks scaling, fitted-model
gradients, held-out probabilities, pairwise AUROC, threshold-block average
precision, tie-adjusted top-K selection and Brier scores. Python hashlib checks
the source, every bundle artifact and executable identity.

Negative tests cover target overlap, group leakage, held-out-label changes,
missing scores, untested outcomes, invalid schemas, malformed CSV, source
tampering, changed imported records even after their recorded hash is updated,
wrong executable identity, extra files, duplicate IDs, nonfinite scores,
single-class training and existing output directories.

The source-fetch helper is tested with mocked response bytes. Its real network
request was not executed successfully in this environment. The tests verify
that an incorrect pinned Git-blob identity or a failed request publishes nothing,
and that existing sources/plans remain unchanged.

## Reproduce

```sh
python3 Tools/BinderBenchmark/check_native.py
python3 Tools/BinderBenchmark/test_tools.py -v
```

[qualification.json](qualification.json) records the exact checked source SHA-256
values, toolchain, counts, campaign dimensions and unexecuted boundaries. This is
an execution record, not a cryptographically signed attestation. Re-run after
source changes. It does not transfer evidence to different implementations.

## Scientific boundary

**No experimentally measured improvement has been established.** The published
CSV was not scored; the synthetic outcomes test software behavior only. This
increment supplies the comparison baseline and traceable evaluator, not Numi
physics features, trained biomolecular model inference or a qualified molecular
preparation/simulation pipeline. Whole-target holdout plus exact-sequence-group
purging does not establish independence from homologous targets or designs.

Next is an explicit published-table run, followed by structural/physical features
from existing molecular owners on the same matched candidates. Freeze the plan
before viewing held-out outcomes, preserve separate assay endpoints and compare
against both the fixed model score and a score-only ensemble.
