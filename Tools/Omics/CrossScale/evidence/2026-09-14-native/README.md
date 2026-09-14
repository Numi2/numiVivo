# Native cross-scale evidence contract qualification (2026-09-14)

The source revision `0ca6dcb8c52d714ccd5ba5bdde81bb791395d2fc` was built and run on the physical Apple M4 Pro Mac mini. `swift test --filter VivoCrossScaleEvidenceTests` passed all three focused tests.

The contract validates the declared research path `variant → regulation → RNA/cell state → protein → molecular mechanism → reaction/kinetics → cellular phenotype → tissue prediction`. It requires adjacent stage order, source and model provenance appropriate to each evidence class, and source/model/validation fingerprints plus held-out observations before a link can be marked qualified. Hypothesis and unavailable links remain explicit and cannot authorize an outcome report.

This is native software and evidence-admission qualification only. The tests use synthetic fingerprints and do not run a public biological dataset, a trained regulatory/protein/phenotype model, a live cross-scale graph, or a tissue/clinical outcome prediction. No current NumiVivo graph has seven independently qualified biological links; the biological outcome claim remains open. `execution.json` binds the command, host, source hashes and test result; `manifest.json` binds every receipt file.
