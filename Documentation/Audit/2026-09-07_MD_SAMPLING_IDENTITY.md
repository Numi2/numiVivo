# MD numerical identity in sampling adapters

The constant-pH and nuclear Metal adapters previously identified their sources,
configuration and optional force provider without identifying the executable MD
numerics. A generic sampler checkpoint could therefore retain an unchanged
Hamiltonian hash across an MD numerical-contract change. Constant-pH propagation
then constructs a fresh MD checkpoint from its physical state; the current MD
checkpoint gate alone cannot distinguish that old outer sampler state.

The Metal-specific constant-pH Hamiltonian hash now includes
`VivoMDExecutionIdentity.current`. The NCMC factory already forwards endpoint
Hamiltonian hashes into its switch-engine hash, and the generic NCMC continuation
validator checks both endpoint and engine identities before any evaluator runs.
Its callback also checks endpoint identity before constructing an MD runtime.
No NCMC algorithm or generic sampler schema changes are needed.

Both nuclear Metal identities include the same contract: the serializable
`VivoNuclearMetalSpecification.definition()` identity, and the direct
`VivoNuclearPotential.metal` evaluator identity. Both are necessary because the
serializable factory deliberately returns its prepared definition around the
direct evaluator. Existing ring-polymer admission and reactive label/model
definition checks then reject stale Metal definitions. Self-contained historical
artifacts can still validate against their historical definitions; this does not
authorize continuation with a different current Metal evaluator.

`MetalSamplingIdentityTests` reconstructs the exact prior identities with the
contract key absent. It proves that these old checkpoints remain structurally
valid for matching generic states but reject current Metal state identities,
isolates the NCMC endpoint and engine checks, rejects old ring admission and force
labels, and checks the direct Metal evaluator hash. Current-identity controls
remain admitted. These are host-only identity and admission tests: no force,
stochastic, GPU or scientific-conformance result is inferred from them.

Generic CPU potential identities and checkpoint schemas are unchanged. Existing
checkpoint data are not relabelled or migrated to a new numerical contract.
