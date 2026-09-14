# Native cross-scale CLI qualification (2026-09-14)

The source revision `85a098c30f4d906a22c4d75a236cd8d98efa2ac5` was built on the physical Apple M4 Pro Mac mini. `swift test --filter VivoCrossScaleEvidenceTests` passed all three focused tests. The `numivivo cross-scale-assess` smoke command emitted the canonical incomplete assessment with exit 0; `--require-qualified` preserved that assessment and returned exit 2 at the first unqualified boundary (`variant → regulation`).

The command reads a bounded `VivoCrossScaleEvidenceGraph` through descriptor-rooted I/O. Existing parent symlinks such as macOS `/tmp` are resolved before the leaf is admitted with `O_NOFOLLOW`; the native run verifies that behavior.

This is native software and evidence-admission qualification only. The synthetic graph contains no biological observations, trained predictor, phenotype, tissue or clinical outcome. It therefore does not establish that NumiVivo can predict biological outcomes.
