# HIRISA integration failure targeting

This artifact turns the existing frozen annotation-retention results into an
explicit target ledger. It does not fit an integration model and does not change
the default global correction.

The complete HIRISA global and condition-stratified results contain the same 471
author-label/stratum comparisons. The ledger only calls a comparison
'supported_control_sensitive' when the pre-existing diagnostic had sufficient
training and query cells in every positive donor fold and a baseline advantage
over erasure. That is a diagnostic support condition, **not** evidence that the
cells are biologically comparable across donors or that correction should be
applied.

On the frozen published inputs:

| Supported, control-sensitive comparisons | Native global failures | Condition-stratified failures | Shared | New under candidate | Resolved under candidate |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 146 | 37 | 40 | 29 | 11 | 8 |

The 325 unsupported comparisons remain in the output with explicit reasons
('insufficient_support', 'insufficient_control_sensitivity', or both). They are
not averaged away or converted into an alignment claim.

The condition candidate uses a separately hashed freeze because it adds its
candidate score matrix to the frozen matrix registry. The copied native and
candidate freeze documents are verified to agree on every other field and on
every pre-existing matrix; the only permitted addition is the pinned
'condition-f1' score matrix. This preserves the original retention design while
making the candidate input explicit.

The next method experiment should start from the 29 shared and 11 newly
introduced supported failures, checking donor-fold behavior and local
composition/state context before proposing a constrained correction. The ledger
does not establish a safe biological alignment domain, does not qualify a
corrected representation, and does not provide phenotype or outcome evidence.

## Reproduce

Use the exact archived score diagnostics, not a newly fitted integration:

\`\`\`sh
python Tools/Omics/Benchmarks/HIRISA/summarize_annotation_retention_failures.py \\
  --ledger Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-11-annotation-retention/study/ledger.json.gz \\
  --native Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-11-annotation-retention/study/native.json.gz \\
  --candidate Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-14-condition-stratified/condition-f1-vs-baseline.json.gz \\
  --native-freeze Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-15-failure-targeting/native-freeze.json \\
  --candidate-freeze Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-15-failure-targeting/condition-f1-freeze.json \\
  --output Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-15-failure-targeting/failure-targets.json
\`\`\`

The generated ledger is machine-readable and captures hashes of all three input
files. See its companion manifest for the exact script and output hashes.
