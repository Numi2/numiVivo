# Digital tissue homeostasis reference

This project exercises the sequence-free computational layers of NumiVivo on Apple silicon. It is a software integration example, not a biological construct, dosing protocol, clinical model, or experimental instruction.

The executable performs one bounded workflow:

1. Defines three abstract physiological compartments and one abstract mediator.
2. Compiles reversible partition coefficients into a mass-conserving Metal model.
3. Advances the model through a transactional RK2 step.
4. Requests before/after amount totals in the step certificate.
5. Packages model and state into a fingerprinted checkpoint envelope.
6. Expands a deterministic Halton parameter campaign with reproducible replicate seeds.
7. Partitions a reaction hypergraph into exact, continuous, and spatial authorities.
8. Evaluates a bounded affine surrogate through the authoritative-fallback gate.
9. Generates the canonical NumanX, NumiTissue, NumiBrain, and NumiVivo channel declarations.
10. Prints one sorted JSON summary.

## Run on Apple silicon

From this directory:

```bash
swift run -c release
```

The root NumiVivo package is referenced through a local Swift Package Manager dependency. A Metal-capable macOS 15 or newer system is required.

No generated checkpoint is written to disk by default. The example constructs it in memory and reports its fingerprint and byte count.

## Expected output structure

```json
{
  "campaignFingerprint": "…",
  "candidateCount": 16,
  "checkpointFingerprint": "…",
  "couplingChannelCount": 30,
  "finalConcentrations": ["…"],
  "hybridModes": ["exactSSA", "deterministicRK2", "spatialSplit"],
  "jobCount": 64,
  "partitionModelFingerprint": "…",
  "partitionStep": { "disposition": "committed" },
  "surrogateAccepted": true,
  "surrogateReason": "accepted"
}
```

Exact numerical values depend on the compiled Metal implementation and device. The model fingerprints, candidate ordering, replicate seeds, campaign ledger, and checkpoint identity are deterministic for identical source artifacts.

## Files

- `Package.swift` — standalone example package.
- `Sources/DigitalTissueHomeostasis/main.swift` — executable integration path.
- `program.json` — declarative molecular-control program.
- `mechanism-problem.json` — abstract mechanism requirements.
- `mechanism-library.json` — sequence-free candidate mechanism records.
- `population-model.json` — spatial phenotype-population source model.
- `partition-model.json` — physiological partition source model.
- `exact-ssa-model.json` — discrete reaction source model.
- `campaign-definition.json` — high-level campaign source.
- `coupling-contract.json` — minimal custom coupling contract.

The JSON source documents are intentionally separate from generated fingerprints and binary packs. Compiler-produced artifacts should be stored in a run bundle rather than edited manually.
