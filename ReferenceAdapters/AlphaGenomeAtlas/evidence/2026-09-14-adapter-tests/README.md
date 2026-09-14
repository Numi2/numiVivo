# AlphaGenome adapter contract verification (2026-09-14)

The pinned adapter sources at commit `bfafd4bbd1d3dba528136a0a7b06440c87475394`
compile with Python 3 and pass all eight offline contract tests. The run uses
no AlphaGenome API key, live provider, pVACtools installation or biological
outcomes.

`execution.json` records the exact source hashes and commands. `test.log`
retains the verbose test output, and `manifest.json` binds both files. The
tests cover request and reference-allele identity, usage and no-training gates,
explicit available/no-data/failed outcomes, duplicate rejection, prepared-job
input hashes, public-reference gating, and a fake RegTools+pVACsplice bundle.

This verifies the adapter boundary only. It does not qualify live Atlas access,
variant-effect accuracy, a regulatory link to RNA/cell state, or any clinical,
phenotypic or treatment outcome.
