# Cross-scale evidence contract

`VivoCrossScaleEvidenceGraph` is the native ownership boundary for the long-term
variant-to-tissue research path:

`variant → regulation → RNA/cell state → protein → molecular mechanism → reaction/kinetics → cellular phenotype → tissue prediction`

The graph is an evidence ledger, not an outcome model. Each node carries a
quantity, unit, context fingerprint and evidence class (`measured`, `predicted`,
`simulated` or `unavailable`). Each link must connect adjacent stages and retain
its source fingerprints, model identity, limitations and validation identity.

The validator enforces these rules:

- stage order is fixed and links cannot skip a boundary or point backwards;
- measured nodes require a source and cannot carry a predictive model;
- predicted and simulated nodes require both source and model fingerprints;
- a `qualified` link requires source, model and validation fingerprints plus at
  least one held-out observation;
- a `hypothesis` link is retained as a hypothesis and cannot authorize an
  outcome report; and
- an `unavailable` link carries no invented evidence and is reported as a gap.

`assess()` reports the seven required boundaries as `incomplete`,
`hypothesisOnly` or `qualifiedResearchPath`. The latter is the only state that
passes `requireQualifiedResearchPath()`, and it means only that the declared
research evidence path is complete. It does not grant clinical, treatment or
patient-specific authorization.

The source contract and focused regression tests live in
`Sources/NumiVivoKit/Omics/VivoCrossScaleEvidence.swift` and
`Tests/NumiVivoIntegrationTests/VivoCrossScaleEvidenceTests.swift`. The tests
use synthetic fingerprints to exercise validation and serialization; they are
software contract evidence. No current NumiVivo graph has seven independently
qualified links, so the cross-scale biological outcome claim remains open.
