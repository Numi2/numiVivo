# Prepared AMBER system to staged Apple MD

This example uses a user-supplied prepared AMBER topology and restart. No biological parameter set, published simulation result or author data is fabricated here. The generated duration/temperature schedule is illustrative and must be chosen for the scientific problem.

```sh
numivivo amber-import-md complex.prmtop --restart complex.rst7 \
  --structure complex.structure.json --system complex.system.json \
  --state complex.initial.json --mapping complex.mapping.json

numivivo md-protocol-template complex.system.json > complex.protocol.json
numivivo md-protocol-validate complex.protocol.json \
  --system complex.system.json --state complex.initial.json
numivivo md-protocol-run complex.protocol.json \
  --system complex.system.json --state complex.initial.json \
  --store ./complex-artifacts > complex.receipt.json
```

The template performs constrained minimization, Maxwell-Boltzmann velocity initialization, NVT, NPT, and NPT production. State and periodic volume transfer between stages are explicit, configuration-bound records. A required minimization gate blocks continuation when it has not converged. Read the blocker report rather than weakening constraints to force a run through.

The MD execution profile currently admits orthogonal cells. An AMBER truncated octahedron remains a valid imported archival structure, but its skew-cell MD is blocked until a closest-lattice-vector backend is installed. No box geometry is silently replaced. Use an appropriately prepared orthogonal source system for this profile.

Trajectory samples are bounded binary chunks stored with their stage identity and per-frame cell. Numerical checkpoints remain distinct from trajectories. A receipt identifies the immutable restart cursor and UUID-named reference. Resuming that cursor creates a new run reference and does not reset velocities or stochastic step counters. See Documentation/Design/MD_PROTOCOL_WORKFLOW.md for commands, output semantics and failure handling.

Implementation status: source only; commands were not executed during this pass. Apple build, Metal validation and numerical qualification remain user-side tasks.
