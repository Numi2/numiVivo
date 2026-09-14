# Native cross-scale contradiction guard qualification (2026-09-14)

The source revision `c8604b8c923df4e88f02e4265467d69898fb0fcc` was built and run on the physical Apple M4 Pro Mac mini. `swift test --filter VivoCrossScaleEvidenceTests` passed all three focused tests.

This correction makes the evidence ledger reject a graph that declares a gap for a boundary already represented by a link. The assessment therefore cannot silently ignore contradictory missing-evidence declarations.

This remains native software and evidence-admission qualification only. The tests use synthetic fingerprints and do not run a public biological dataset, a trained cross-scale predictor, a live cross-scale graph, or a tissue/clinical outcome prediction. The biological outcome claim remains open.
