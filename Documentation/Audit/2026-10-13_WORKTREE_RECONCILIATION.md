# Recovered workflow proposals and current native verification

Current numerical main `bafff6dd164c449eb6cae4b92f0131dfa438c539` already
contains the archive, conditional encounter, prepared-reaction route, NPT,
PME/thermal qualification and sampling export/CLI proposals, with later ownership
and bounded-reader repairs. The original source workspaces are preserved under
recovery refs and exact staged/unstaged/untracked/unique-ignored archives.

The real CLI checker found that later global help omitted
`molecule-sampling-export-verify`, despite its command and molecule-specific help
remaining available. This change restores discovery and leaves numerical owners
unchanged. The failing check and its source fixture are retained as evidence.

Fresh native validation on Apple M4 passes the 11-test prepared-sampling → fresh
reaction → conditional-encounter route, 41 storage/archive/export host tests,
classical PME (2 tests), multipole PME (3), NPT transition (1), thermal (2 tests, 3 argument cases), and
accepted sampling prefix/resume (1). The prepared reaction test takes about
2,481 seconds; this is a completed fresh electronic calculation. Current real
sampling export/verify/connected-workflow checks pass all 87 checks, and the
workflow CLI campaign passes 61 checks after materializing its declared inputs.
The receipt records original/repaired binary and log SHA-256 identities.

These are the scopes stated by the existing tests. Budget-exhausted sampling is
not a converged ensemble; synthetic harmonic routing is not physical H2;
force-free thermal and ideal volume transitions are not interacting NPT
qualification; finite basis/reaction fixtures are not biological accuracy.
Historical failed and contradicted attempts remain available. No performance or
cross-hardware improvement follows from this verification.
