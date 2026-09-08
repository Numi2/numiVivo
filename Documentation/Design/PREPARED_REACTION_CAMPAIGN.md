# Verified geometry seeds and conditional encounter observables

These interfaces connect preparation, accepted molecular sampling, fresh reaction qualification and a time-dependent observable. They reuse the existing preparation, Metal MD, snapshot mapper, electronic/path, workflow and kinetic authorities.

## Accepted state authority

`VivoMolecularSamplingExporter.verifiedExport` performs the same bounded archive and immutable receipt verification as `verify`, returning a non-Codable `VivoVerifiedMolecularSamplingExport`. Only the exporter can construct it. It retains the exact receipt, verification summary, accepted checkpoint bytes and rooted store. Decoding the serializable verification summary does not grant this authority.

## Geometry transfer

`VivoSamplingReactionSeedRequest` binds a source export, a destination reaction request, explicit atom mappings and a bounded model-transfer statement. Destinations are fresh `.qualify` or `.connectedReaction` requests. Each assignment names either the standalone qualification, connected saddle, or an endpoint identifier/component index, and maps every destination nucleus to a distinct source chemical atom.

The assembler uses the existing atom-to-particle snapshot mapping and converts nanometres to Bohr. It changes only initial nuclear coordinates. Atomic numbers must agree; masses, basis, electron counts, solver, optimizer, thermochemical temperature and standard state remain the supplied destination settings. Periodic source geometry and assigned external point-charge environments are rejected because this interface defines no unwrapping, partitioning or environmental transformation. Already qualified points and Hessians cannot be destinations.

The publication retains source and destination hashes, mapped particles, masses, temperature/model contexts, accepted-step/status evidence and the transfer statement. Its typed `request` output is a normal `vivo.reaction-calculation-request` workflow artifact. Fresh numerical qualification remains mandatory. Read-only verification reopens the source export and reconstructs the exact immutable assembly receipt without consulting mutable caches or regenerating missing artifacts.

## Qualified rate to observable

`VivoConditionalEncounterRequest` identifies one tagged component of the declared reactant endpoint and assigns every other mapped component to an externally maintained reservoir. Each component binding includes its index, ordered global atom identities and canonical qualified-point hash. Reservoir values are explicit partial pressures in Pa or ideal concentrations in mol/L, consistent with the rate's standard state. Temperature must match the source. Missing, duplicate, reidentified or dimensionally incompatible components reject.

The hazard is `k × product(reservoir values)`; it has units s⁻¹. A unimolecular source needs no reservoir. Zero reservoir values yield zero hazard, while a logarithmic derivative with respect to a zero reservoir is undefined and rejects. Products are calculated in logarithmic form; unrepresentable positive hazards reject explicitly.

The existing `VivoLinearKineticSensitivity` evolves the one-state absorbing generator `[-hazard]` from initial probability `[1]`. It supplies survival, first-event probability, flux, hazard and local natural-log-rate/optional reservoir sensitivities. Its existing work limits and truncation accounting remain authoritative. A condition-building preflight can check these bounds, but cannot qualify a decoded reaction.

`VivoConditionalEncounter.calculate` freshly validates the complete source TST result through the existing connectivity, thermochemistry and rate reconstruction before returning success, including for zero reservoirs. The workflow operation consumes the exact `vivo.reaction-calculation-result` input artifact plus the explicit encounter request, and retains the source payload hash and original transmission declaration. Output validation reconstructs the result through the same authorities.

These are conditional first-event predictions. No finite-reservoir depletion, reverse kinetics, occupancy, equilibrium-population, efficacy, independent transmission or physical uncertainty claim is added. A partial classical sampling prefix remains a geometry seed; it is not promoted to a QM ensemble or averaged reaction rate.

## Executable example

The [published preparation-to-observable example](../../Examples/prepared-reaction/README.md) supplies every model and bath input. The public runner retains separate records for each stage and checks failure cases as well as the successful route. Native and CLI evidence must identify the source, executable and actual execution scope before this example is described as verified.
