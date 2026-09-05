# Staged Apple-native MD

This example connects the existing Wave B MD engine to a single reproducible
minimization → NVT → NPT → production schedule. It does not bundle an AMBER
parameter set or a prepared biological system. Supply a compatible prepared
`system.prmtop` and `system.rst7` pair that you are entitled to use.

The supplied base configuration uses a 1 fs integration step, 300 K target,
PME electrostatics, 1 nm real-space cutoff and a 0.15 nm neighbor skin. Those are
explicit example inputs, not a recommendation or qualification for every model.
The periodic cell and constrained model must satisfy the backend's capability
checks. Step counts below demonstrate the interface; they do not certify
thermodynamic equilibration or convergence.

After building the CLI on the user's Apple-silicon machine, from the repository root:

```sh
numivivo amber-import-md system.prmtop --restart system.rst7 \
  --structure prepared.structure.json --system prepared.system.json \
  --state prepared.initial.json --mapping prepared.mapping.json

numivivo md-protocol-template \
  --config Examples/md-staged-protocol/base-md.json \
  --nvt-steps 50000 --npt-steps 50000 --production-steps 500000 \
  --sample-every 1000 --observe-every 1000 --checkpoint-every 1000 \
  --output protocol.json

numivivo md-protocol-check prepared.system.json \
  --state prepared.initial.json --protocol protocol.json \
  --output protocol.capabilities.json

numivivo md-protocol-run prepared.system.json \
  --state prepared.initial.json --protocol protocol.json \
  --checkpoint protocol.checkpoint.json \
  --trajectory trajectory-001.nvmd --output protocol.result.json

numivivo md-trajectory-inspect trajectory-001.nvmd \
  --output trajectory-001.inspection.json
```

The generated protocol initializes constrained thermal velocities once after
minimization. It does not reinitialize them at every stage. Subsequent stages
preserve positions, velocities, box and physical time while using explicit new
configuration identities and distinct thermostat seeds.

An interrupted run can continue from its last published protocol checkpoint:

```sh
numivivo md-protocol-run prepared.system.json \
  --restore protocol.checkpoint.json --protocol protocol.json \
  --checkpoint protocol.continued.checkpoint.json \
  --trajectory trajectory-002.nvmd --output protocol.continued.result.json
```

Use a new trajectory segment. The archive is finalized atomically and is not
opened for in-place append during resume. The restarted segment may repeat the
accepted boundary frame; its metadata identifies the stage, local step and
physical time. It does not repeat the completed dynamics or the thermalization.

A protocol checkpoint is not interchangeable with a plain `md-run` checkpoint.
It also binds the stage definition and progress. Keep the original protocol
unchanged for resume. Starting a new experiment with changed settings is a
separate operation, not a way to overwrite checkpoint identity.

Coordinate archives store packed FP32 triples with SHA-256 frame checksums and
a final record-chain footer. `--include-velocities true` also stores velocities.
Frames stream incrementally rather than accumulating every position in RAM.
Observations are collected independently of the coordinate sampling interval.
The JSON result contains stage summaries, observations, checkpoint identity and
archive receipt, not a second complete copy of the coordinate trajectory.

Preflight rejects output/input path collisions and existing outputs unless
`--force` is explicit. Interval checkpoint replacement is atomic. On a failed
required minimization or rejected MD step, the protocol stops and writes its
latest available accepted state rather than claiming all stages completed.

## Qualification still required

No commands above were executed during this source-development continuation.
Apple compilation, Metal validation, native force comparison, PME error analysis,
constraint/minimization behavior, ensemble sampling and checkpoint equivalence
remain separate gates. Read
`Documentation/Design/WAVE_B_PROTOCOLS_AND_TRAJECTORIES.md` for precise format,
transaction boundaries, and inherited numerical limitations.
