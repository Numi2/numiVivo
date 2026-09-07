import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct PlatformSnapshotMappingHelperTests {
    private struct Fixture {
        let source: VivoMolecularStructureDocument
        var system: VivoClassicalSystem
        var checkpoint: VivoMDCheckpoint
    }

    private func fixture() throws -> Fixture {
        let recipe = try VivoWorkflowTemplates.molecularDynamics()
        func inline<T: Decodable>(_ name: String, _ type: T.Type) throws -> T {
            guard let artifact = recipe.artifacts.first(where: { $0.identifier == name }),
                  case .json(_, let payload) = artifact.source else {
                throw VivoChemistryError.invalid("missing snapshot helper fixture input")
            }
            return try VivoPlatformOperations.decode(type, payload)
        }
        let source = try inline("structure", VivoMolecularStructureDocument.self)
        var system = try inline("system", VivoClassicalSystem.self)
        system.particles[0].atomIndex = 1
        system.particles[1].atomIndex = 0
        system.particles.append(.init(index: 2, atomIndex: nil, typeIdentifier: "snapshot-site",
            role: .virtualSite, massDa: 0, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0))
        system.linearVirtualSites = [.init(siteParticle: 2, parentParticles: [0, 1], weights: [0.5, 0.5])]
        let checkpoint = VivoMDCheckpoint(systemFingerprint: try system.fingerprint(),
            configurationFingerprint: try VivoCanonicalJSON.fingerprint(Data("snapshot-helper-fixture".utf8)),
            acceptedStep: 37, timePS: 0.125,
            positionsNM: [.init(4.5, -0.25, 0), .init(1.25, 2, 3), .init(2.875, 0.875, 1.5)],
            velocitiesNMPerPS: [.init(0.5, -0.25, 0), .init(-0.5, 0.25, 1), .init(0, 0, 0.5)],
            periodicCell: .init(a: .init(4, 0, 0), b: .init(0.5, 4, 0), c: .init(0, 0.5, 4)))
        return .init(source: source, system: system, checkpoint: checkpoint)
    }

    @Test func mapsOnlyChemicalAtomsAndPreservesAcceptedGeometryAndLineage() throws {
        let f = try fixture()
        let snapshot = try VivoPlatformSnapshotOperations.map(source: f.source, system: f.system, checkpoint: f.checkpoint)
        let checkpointID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(f.checkpoint))
        #expect(snapshot.mapping.atomToParticle == [1, 0])
        #expect(snapshot.frame.positionsNM == [f.checkpoint.positionsNM[1], f.checkpoint.positionsNM[0]])
        #expect(snapshot.frame.velocitiesNMPerPS == [f.checkpoint.velocitiesNMPerPS[1], f.checkpoint.velocitiesNMPerPS[0]])
        #expect(snapshot.frame.step == 37 && snapshot.frame.timePS == 0.125)
        #expect(snapshot.frame.periodicCell == f.checkpoint.periodicCell)
        #expect(snapshot.document.structure.periodicCell == f.checkpoint.periodicCell)
        #expect(snapshot.document.structure.atoms == f.source.structure.atoms)
        #expect(snapshot.document.structure.bonds == f.source.structure.bonds)
        #expect(snapshot.document.structure.conformers.count == 1)
        #expect(snapshot.document.structure.conformers[0].identifier == "md-step-37")
        #expect(snapshot.document.structure.conformers[0].positionsNM == snapshot.frame.positionsNM)
        #expect(snapshot.document.sourceFingerprint == checkpointID)
        #expect(snapshot.mapping.checkpointPayload == checkpointID)
        #expect(snapshot.mapping.sourceStructure == f.source.structureFingerprint)
        #expect(snapshot.mapping.classicalSystem == f.checkpoint.systemFingerprint)
        #expect(snapshot.mapping.snapshotStructure == snapshot.document.structureFingerprint)
        #expect(snapshot.mapping.step == snapshot.frame.step && snapshot.mapping.timePS == snapshot.frame.timePS)
        #expect(snapshot.document.structure.metadata["numivivo.snapshot.sourceStructure"] == f.source.structureFingerprint.hex)
        #expect(snapshot.document.structure.metadata["numivivo.snapshot.classicalSystem"] == f.checkpoint.systemFingerprint.hex)
        #expect(snapshot.document.structure.metadata["numivivo.snapshot.checkpointPayload"] == checkpointID.hex)
        #expect(snapshot.mapping.convention == "atom order preserved; non-atomic particles excluded; coordinates nm, velocities nm/ps; original identity and checkpoint retained; no unwrapping, missing-atom repair, electron assignment or QM/MM approximation inferred")
    }

    @Test func workflowKeepsCanonicalPayloadIdentityAndIdenticalOutputBytes() async throws {
        let f = try fixture()
        let canonicalCheckpoint = try VivoCanonicalJSON.encode(f.checkpoint)
        let originalCheckpoint = Data(" \n".utf8) + canonicalCheckpoint + Data("\n ".utf8)
        let canonicalID = try VivoCanonicalJSON.fingerprint(canonicalCheckpoint)
        #expect(try VivoCanonicalJSON.fingerprint(originalCheckpoint) != canonicalID)
        #expect(try VivoCanonicalJSON.decode(VivoMDCheckpoint.self, from: originalCheckpoint) == f.checkpoint)
        let snapshot = try VivoPlatformSnapshotOperations.map(source: f.source, system: f.system,
            checkpoint: f.checkpoint, budget: .init())
        let expected = ["structure": try snapshot.document.canonicalData(),
                        "frame": try VivoCanonicalJSON.encode(snapshot.frame),
                        "mapping": try VivoCanonicalJSON.encode(snapshot.mapping)]
        let operation = VivoPlatformSnapshotOperations.definition(implementationFingerprint: canonicalID).operation
        var inputs = ["structure": try f.source.canonicalData(), "system": try VivoCanonicalJSON.encode(f.system),
                      "checkpoint": canonicalCheckpoint]
        let canonicalOutputs = try await operation.execute(.object([:]), inputs, .init())
        #expect(canonicalOutputs == expected)
        inputs["checkpoint"] = originalCheckpoint
        let originalOutputs = try await operation.execute(.object([:]), inputs, .init())
        #expect(originalOutputs == canonicalOutputs)
        try operation.validateOutputs(.object([:]), inputs, originalOutputs, .init())
        let mapping = try VivoPlatformOperations.input(VivoWorkflowSnapshotMapping.self, "mapping", originalOutputs)
        #expect(mapping.checkpointPayload == canonicalID)
    }

    @Test(arguments: ["duplicate", "missing", "virtual"])
    func rejectsInvalidChemicalAtomOwnership(_ defect: String) throws {
        var f = try fixture()
        switch defect {
        case "duplicate": f.system.particles[1].atomIndex = 1
        case "missing": f.system.particles[0].atomIndex = nil
        default: f.system.particles[2].atomIndex = 0
        }
        f.checkpoint.systemFingerprint = try f.system.fingerprint()
        do {
            _ = try VivoPlatformSnapshotOperations.map(source: f.source, system: f.system, checkpoint: f.checkpoint)
            Issue.record("invalid snapshot chemical atom ownership was accepted")
        } catch VivoChemistryError.invalid(let message) {
            #expect(message.contains("snapshot"))
        }
    }

    @Test(arguments: [false, true])
    func rejectsForeignCheckpointOrSource(_ foreignSource: Bool) throws {
        var f = try fixture()
        if foreignSource {
            f.system.structureFingerprint = f.checkpoint.configurationFingerprint
            f.checkpoint.systemFingerprint = try f.system.fingerprint()
        } else {
            f.checkpoint.systemFingerprint = f.checkpoint.configurationFingerprint
        }
        do {
            _ = try VivoPlatformSnapshotOperations.map(source: f.source, system: f.system, checkpoint: f.checkpoint)
            Issue.record("foreign snapshot state was accepted")
        } catch VivoChemistryError.invalid(let message) {
            #expect(message == "snapshot checkpoint, classical system and source structure differ")
        }
    }

    @Test func preservesCheckpointAndResourceAdmission() throws {
        var f = try fixture()
        do {
            _ = try VivoPlatformSnapshotOperations.map(source: f.source, system: f.system,
                checkpoint: f.checkpoint, budget: .init(maximumBytes: 1))
            Issue.record("snapshot helper ignored its resource budget")
        } catch VivoChemistryError.resourceLimit {}
        f.checkpoint.numericalContract = nil
        do {
            _ = try VivoPlatformSnapshotOperations.map(source: f.source, system: f.system, checkpoint: f.checkpoint)
            Issue.record("snapshot helper accepted a checkpoint without its numerical contract")
        } catch VivoArtifactValidationError.incompatible {}
        f.checkpoint.numericalContract = VivoMDExecutionIdentity.current
        f.checkpoint.positionsNM[0].x = 0.1
        do {
            _ = try VivoPlatformSnapshotOperations.map(source: f.source, system: f.system, checkpoint: f.checkpoint)
            Issue.record("snapshot helper accepted state requiring another FP32 rounding")
        } catch VivoArtifactValidationError.invalid {}
    }
}
