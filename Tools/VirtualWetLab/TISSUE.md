# Shared tissue representation

`VivoTissueSpecimen` is owned by NumiVivo, independent of NumiLab's renderer.
It retains a hierarchy of molecules, cells, neighborhoods, regions and specimen
entities; named biological units; source hashes; original coordinate systems;
and sparse measurements with independent modality, units, time and evidence.
RNA, protein, morphology, extracellular fields and image references fit the same
contract. A missing sparse entry means zero only when explicitly declared.

Evidence states distinguish MEASURED, SIMULATED — VALIDATED DOMAIN, SIMULATED —
OUT OF DISTRIBUTION, MODEL INFERENCE, HYPOTHESIS and UNAVAILABLE. Model output
requires model identity; validated-domain output requires an explicit domain.
Uncertainty and source identity are mandatory. Passing schema checks does not
validate a biological model or independently authenticate a source's claims.

Build the bounded native validator with `build_tissue.py --output NEW_DIRECTORY`,
then use `numivivo-tissue specimen.json`. `check_tissue.py --binary EXECUTABLE
--output NEW_DIRECTORY` exercises real native admission and rejection. The root
Swift package also exposes the executable. No existing Omics or transport model
has been replaced.

`VivoMechanobiologyParticipant` defines the next cross-owner boundary: private
prepare, mutually accepted release and explicit abort/checkpoint verification.
Channels carry quantity, units, entity IDs, model identity and evidence. It does
not execute bidirectional Matter coupling. Implementations must retain the native
prepare/release authority illustrated by `VivoMolecularPhysiologyCoordinator`;
a failed distributed release must quarantine the experiment, not claim rollback.
The current Matter-to-transport adapter remains a one-way frozen geometry bridge.

The exchange channel now carries named vector components and per-entity values,
units, model/source identity and evidence. Proposals reject unsupported channel
directions, unavailable values, malformed shapes and invalid clocks. Prepared-pair
admission checks both owner checkpoints, transaction identity, identical end time
and declared convergence tolerance before a caller may request release. This
structural gate cannot prove native convergence, atomic commit or rollback.
`check_mechanobiology.py --output NEW_DIRECTORY` compiles an external Swift client
and exercises valid admission plus seven rejection cases. Public initializers
allow both NumiLab and other native owners to construct shared tissue values.

The real spatial experiment retains its original v1 validator binary and receipt
for exact replay. The v2 interface build is separate; it does not replace that
sealed runtime or alter the retained biological comparison.
