import Foundation

public struct VivoWorkflowSnapshotMapping: Codable, Sendable, Equatable {
    public let sourceStructure: VivoFingerprint
    public let classicalSystem: VivoFingerprint
    /// Digest of the checkpoint's canonical payload, not its output envelope.
    public let checkpointPayload: VivoFingerprint
    public let snapshotStructure: VivoFingerprint
    public let atomToParticle: [UInt32]
    public let step: UInt64
    public let timePS: Double
    public let convention: String
}

/// Accepted particle state to the existing atom-based structure/trajectory IR.
/// Virtual particles never become extra chemical atoms. A changed conformer gets
/// a new structure identity; the original topology/system identities remain in
/// the mapping artifact. Periodicity is preserved, not silently unwrapped.
public enum VivoPlatformSnapshotOperations {
    public static func definition(implementationFingerprint id: VivoFingerprint) -> VivoWorkflowDefinition {
        VivoPlatformOperations.pure(identifier: "vivo.platform.md-snapshot", id: id,
            inputs: ["structure": "vivo.molecular-structure-document", "system": "vivo.classical-system", "checkpoint": "vivo.md-checkpoint"],
            outputs: [.init(name: "structure", kind: "vivo.molecular-structure-document"),
                      .init(name: "frame", kind: "vivo.trajectory-frame"),
                      .init(name: "mapping", kind: "vivo.md-snapshot-mapping")],
            summary: "Accepted MD checkpoint to atom-mapped structure and trajectory frame; retains periodicity and source lineage.",
            configure: VivoPlatformOperations.empty, calculate: { _, inputs, budget in
                let source = try VivoPlatformOperations.input(VivoMolecularStructureDocument.self, "structure", inputs)
                let system = try VivoPlatformOperations.input(VivoClassicalSystem.self, "system", inputs)
                let checkpoint = try VivoPlatformOperations.input(VivoMDCheckpoint.self, "checkpoint", inputs)
                try source.validate(); try budget.validate()
                let n = source.structure.atoms.count
                guard n > 0, n <= Int(UInt32.max) else { throw VivoChemistryError.invalid("MD snapshot atom count") }
                _ = try budget.elements([n, 3], simultaneousArrays: 16)
                try VivoClassicalSystemValidator.validate(system, atomCount: UInt32(n))
                try checkpoint.validate(particleCount: system.particles.count)
                let systemID = try system.fingerprint()
                guard source.structureFingerprint == system.structureFingerprint,
                      checkpoint.systemFingerprint == systemID else {
                    throw VivoChemistryError.invalid("snapshot checkpoint, classical system and source structure differ")
                }
                var mapping = [Int](repeating: -1, count: n)
                for particle in system.particles {
                    if particle.role == .atom {
                        guard let atom = particle.atomIndex, atom < UInt32(n), mapping[Int(atom)] == -1 else {
                            throw VivoChemistryError.invalid("snapshot physical particle/atom mapping is absent or duplicated")
                        }
                        mapping[Int(atom)] = Int(particle.index)
                    } else if particle.atomIndex != nil {
                        throw VivoChemistryError.invalid("non-atomic particle cannot own a snapshot chemical atom")
                    }
                }
                guard mapping.allSatisfy({ $0 >= 0 }) else { throw VivoChemistryError.invalid("snapshot has unmapped chemical atoms") }
                let frame = VivoTrajectoryFrame(step: checkpoint.acceptedStep, timePS: checkpoint.timePS,
                    positionsNM: mapping.map { checkpoint.positionsNM[$0] },
                    velocitiesNMPerPS: mapping.map { checkpoint.velocitiesNMPerPS[$0] }, periodicCell: checkpoint.periodicCell)
                try frame.validate(atomCount: n)
                let checkpointID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(checkpoint))
                var structure = source.structure
                structure.conformers = [.init(identifier: "md-step-\(checkpoint.acceptedStep)", positionsNM: frame.positionsNM)]
                structure.periodicCell = frame.periodicCell
                structure.metadata["numivivo.snapshot.sourceStructure"] = source.structureFingerprint.hex
                structure.metadata["numivivo.snapshot.classicalSystem"] = systemID.hex
                structure.metadata["numivivo.snapshot.checkpointPayload"] = checkpointID.hex
                let document = try VivoMolecularStructureDocument(structure: structure, sourceFingerprint: checkpointID)
                let lineage = VivoWorkflowSnapshotMapping(sourceStructure: source.structureFingerprint,
                    classicalSystem: systemID, checkpointPayload: checkpointID, snapshotStructure: document.structureFingerprint,
                    atomToParticle: mapping.map(UInt32.init), step: frame.step, timePS: frame.timePS,
                    convention: "atom order preserved; non-atomic particles excluded; coordinates nm, velocities nm/ps; original identity and checkpoint retained; no unwrapping, missing-atom repair, electron assignment or QM/MM approximation inferred")
                return ["structure": try document.canonicalData(), "frame": try VivoCanonicalJSON.encode(frame),
                        "mapping": try VivoCanonicalJSON.encode(lineage)]
            })
    }
}
