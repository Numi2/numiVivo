# Accepted sampling states in research workflows

`VivoMolecularSamplingArchiveReader` selects an exact accepted replica from a completed cross-replica block. `VivoMolecularSamplingExporter` publishes that selection as typed chemistry workflow outputs. Downstream operations receive the original checkpoint, model and explicit atom mapping through the existing artifact store and workflow envelopes.

This interface does not advance MD or turn a partially sampled state into an equilibrated ensemble. A valid accepted prefix can be useful for inspection and further research even when its declared stopping criteria are unmet. The export records that distinction.

## Select, export and verify

```sh
numivivo molecule-sampling-export --store STORE \
  --checkpoint SAMPLING_CHECKPOINT_SHA256 --replica 0 \
  --receipt SAMPLING_RECEIPT_SHA256 --output export.json \
  --checkpoint-output checkpoint.json --snapshot-output snapshot.json \
  --mapping-output mapping.json

numivivo molecule-sampling-export-verify export.json --store STORE \
  --verify-all-payloads --output verification.json
```

`--replica` is an explicit zero-based index. `--checkpoint` names the aggregate sampling cursor, not an individual MD checkpoint. Alternatively, `--reference NAME` resolves one named cursor through the captured store. The runner's reference is `molecular-sampling-<request-sha256>-checkpoint`. Its former slash-containing name was incompatible with the store's reference contract; the runner now generates a valid flat name. Immutable cursor hashes remain the direct selection interface.

`--receipt` is optional and names a stored `molecular-sampling-receipt`. It must bind the selected cursor, request, manifests and diagnostics. Without one, termination is `notRecorded`; a cursor cannot distinguish an active run from a crash, rejection or cancellation. With one, `sourceTermination` retains `converged`, `budgetExhausted`, `rejected` or `cancelled` independently of `declaredObservableCriteriaSatisfied`.

A valid partial export exits 0. `--require-converged` makes satisfaction of the declared consecutive observable criteria a caller policy: unmet criteria exit 65. Invalid inputs also exit 65; cooperative task cancellation exits 130. Process signal behavior is unchanged. Block zero has no accepted production state and cannot be selected. `molecule-sampling-run` retains its existing status convention, including exit 75 for an unconverged terminal result.

Every file destination is checked together before publication, must be outside artifact storage, and must not alias an input or another output. After these collective checks, output directory authority is captured before asynchronous export. Later path replacement cannot redirect these prepared writes. Outputs are no-clobber. Optional files are written individually, followed by the export receipt. A later failure can leave completed optional files or immutable intermediate artifacts; the files are not one filesystem transaction.

## What is verified

Opening a cursor checks all replica checkpoint/manifest/scalar metadata, request identity, replica seeds and configurations, original source clocks, accepted-step and sampling schedules, current numerical contracts. Accepted NPT cell state is authoritative. Selection checks current execution capability of the selected accepted state and additionally verifies the chosen trajectory's complete index chain and its final coordinate payload against the exact checkpoint. `--verify-all-payloads` checks every selected coordinate payload; it does not imply that every other replica's full coordinate archive was inspected.

The reader reconstructs the current diagnostic and the final consecutive-pass suffix from the bound scalar prefixes. It does not remeasure all scalar observables from coordinates or prove that no much earlier stopping opportunity existed. Minimized-state references are retained without revalidating the optimizer history. Declared-observable convergence is not proof of conformational exploration, force-field accuracy, model transfer to QM, or a reaction rate.

Runner continuation shares this validation owner and checks every replica archive with one cumulative budget. It restores captured accepted checkpoints through the existing MD runtime and resumes the verified trajectory writer. It does not rethermalize resumed replicas or reconstruct velocities from coordinate-only frames.

Fresh export verification reopens the source archive, reconstructs the selected outputs and checks the exact immutable workflow receipt and output envelopes. It never executes an operation to regenerate missing outputs. Mutable cache references are irrelevant to this verification. A stronger payload check can verify an export originally made with final-payload validation; the verification reports both recorded and effective scope without changing the original receipt.

## Typed outputs and identity

| Output name | Workflow kind | Meaning |
| --- | --- | --- |
| `checkpoint` | `vivo.md-checkpoint` | Original stored checkpoint bytes, including their original JSON formatting. |
| `source-structure` | `vivo.molecular-structure-document` | Original chemical structure and its identity. |
| `system` | `vivo.classical-system` | Bound classical model and particle ownership. |
| `configuration` | `vivo.md-configuration` | Selected replica configuration, including its seed. |
| `snapshot` | `vivo.md-state-snapshot` | Accepted all-particle coordinates, velocities, cell, step and clock. No new energy or temperature measurement is invented. |
| `structure` | `vivo.molecular-structure-document` | Accepted coordinates in chemical atom order. |
| `frame` | `vivo.trajectory-frame` | Chemical-atom coordinates and velocities with the accepted clock and cell. |
| `mapping` | `vivo.md-snapshot-mapping` | Explicit atom-to-particle map and source/checkpoint/derived identities. |
| `provenance` | `vivo.molecular-sampling-selected-state` | Selection lineage, scope, diagnostics and termination. |

These outputs use `chemistry-output` envelopes. The raw `md-checkpoint` object is preserved with its original kind. The export receipt separately binds task identity, immutable workflow receipt, output-envelope hashes and decoded payload hashes.

The CLI uses the same executable/operating-system implementation fingerprint as general workflow commands. Fresh verification requires that matching implementation. Retain the producing executable and resources for historical verification; another implementation can make a new export from an admitted source prefix.

The snapshot mapper's existing `checkpointPayload` identity hashes the canonical encoding of the decoded state. For valid noncanonical input bytes that hash can differ from the original checkpoint artifact hash. Export provenance retains both as `canonicalCheckpointFingerprint` and `mdCheckpointFingerprint`; neither is substituted for the other.

The mapper preserves periodic cells, excludes non-atomic particles from the chemical frame and rejects ambiguous ownership. It does not unwrap geometry, repair missing atoms, infer electrons or select a basis. Existing `vivo.platform.md-snapshot` and structure-to-electronic operations consume the typed outputs. Periodic electronic work still requires the explicit supported QM/MM interface and partition policy. A selected geometry alone does not supply a qualified saddle, free-energy difference or kinetic model.

## Admission limits

`--read-limits FILE` accepts a bounded JSON document with schema `numivivo.org/molecular-sampling-read-limits/v1`. Optional fields override the library defaults; unknown fields reject. `numivivo molecule-help` lists every field and default. The same option applies to `molecule-sampling-run` only with `--resume`.

```json
{"schema":"numivivo.org/molecular-sampling-read-limits/v1","maximumReadBytes":8589934592,"maximumTrajectoryChunks":200000}
```

Individual request/cursor/checkpoint/diagnostic caps apply before allocation. Cumulative reader admission includes conservative index and coordinate-wire ceilings, with separate scalar, replica, particle, accepted-step and diagnostic-work limits. It is not measured RSS. A null chunk ceiling leaves no separate chunk-count cap; zero rejects nonempty trajectories. Store descriptor limits remain separate.

Named-reference resolution performs an additional integrity read bounded by `maximumCursorBytes`. Strengthening a final-payload export reads one extra manifest and tail link, each capped at 64 KiB, to reconstruct the recorded byte count. Immutable workflow verification has separate checked JSON/base64 wire limits derived from its task contract. Normal workflow producer input semantics are unchanged.

`--maximum-export-bytes BYTES` controls the independent chemistry publication budget, default 268435456 bytes. The budget also bounds the publication task's aggregate input/output payloads. Raising reader limits does not automatically raise this publication budget. None of these limits changes the simulation model.

## Native evidence boundary

`MolecularSamplingNativeBridgeTests` uses a small, explicitly synthetic two-particle harmonic model, two Langevin seeds and deliberately unmet retention criteria. Its success receipt requires real accepted Metal sampling, selection of an actually published earlier cursor, exact continuation to the same final cursor/diagnostics/manifests, terminal resume without advancement, and exact selected state identities. `NUMIVIVO_TEST_ARTIFACTS` retains its store and `NUMIVIVO_TEST_SOURCE_COMMIT` records the tested source.

The process checker consumes that retained native prefix, exercises the public export/verification commands and routes the exported state into explicitly configured electronic operations. Host integrity fixtures separately cover corruption, original byte identity, mapping, limits, cache eviction and cancellation. These are continuation and workflow tests; they do not establish equilibrium statistics, physical H2 parameter accuracy, protein-scale kinetics or performance. Executed results belong in the dated audit and validation record for the tested revision.
