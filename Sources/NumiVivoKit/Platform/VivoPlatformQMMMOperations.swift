import Foundation

public struct VivoWorkflowMDQMMMConfiguration: Codable, Sendable, Equatable {
    public var request: VivoQMMMRegionRequest
    public var solventPolicy: VivoQMMMSolventPromotionPolicy?
    public var cluster: VivoQMMMClusterConfiguration
    public init(request: VivoQMMMRegionRequest, solventPolicy: VivoQMMMSolventPromotionPolicy? = nil,
                cluster: VivoQMMMClusterConfiguration = .init()) {
        self.request = request; self.solventPolicy = solventPolicy; self.cluster = cluster
    }
    public func validate() throws {
        try cluster.validate(); try solventPolicy?.validate()
        guard !request.qmAtomIndices.isEmpty, Set(request.qmAtomIndices).count == request.qmAtomIndices.count,
              (0...100000).contains(request.alphaElectrons), (0...100000).contains(request.betaElectrons),
              request.hydrogenLinkDistancesNM.values.allSatisfy({ $0.isFinite && $0 > 0 && $0 < 1 }) else {
            throw VivoChemistryError.invalid("MD/QM/MM workflow base region")
        }
    }
}

/// Registered beside the isolated snapshot adapter. Neither adapter silently
/// replaces the other's Hamiltonian, periodicity, electron sector or atom mapping.
public enum VivoPlatformQMMMOperations {
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        [VivoPlatformOperations.pure(identifier: "vivo.platform.md-qmmm", id: id,
            inputs: ["structure": "vivo.molecular-structure-document", "system": "vivo.classical-system",
                     "checkpoint": "vivo.md-checkpoint"],
            outputs: [.init(name: "frame", kind: "vivo.md-qmmm-prepared-frame"),
                      .init(name: "system", kind: "vivo.electronic-system"),
                      .init(name: "mapping", kind: "vivo.qmmm-finite-cluster")],
            summary: "Accepted MD checkpoint to molecule-preserving finite QM/MM geometry, solvent promotion and source mapping.",
            configure: { try VivoPlatformOperations.decode(VivoWorkflowMDQMMMConfiguration.self, $0).validate() },
            calculate: { cfg, inputs, budget in
                let configuration = try VivoPlatformOperations.decode(VivoWorkflowMDQMMMConfiguration.self, cfg)
                let document = try VivoPlatformOperations.input(VivoMolecularStructureDocument.self, "structure", inputs)
                let system = try VivoPlatformOperations.input(VivoClassicalSystem.self, "system", inputs)
                let checkpoint = try VivoPlatformOperations.input(VivoMDCheckpoint.self, "checkpoint", inputs)
                let frame = try VivoQMMMCompiler.prepareMD(document: document, system: system, checkpoint: checkpoint,
                    request: configuration.request, solventPolicy: configuration.solventPolicy,
                    configuration: configuration.cluster, budget: budget)
                return ["frame": try VivoCanonicalJSON.encode(frame),
                        "system": try VivoCanonicalJSON.encode(frame.region.electronicSystem),
                        "mapping": try VivoCanonicalJSON.encode(frame.cluster)]
            }),
         VivoPlatformOperations.pure(identifier: "vivo.platform.qmmm-forces", id: id,
            inputs: ["structure": "vivo.molecular-structure-document", "system": "vivo.classical-system",
                     "frame": "vivo.md-qmmm-prepared-frame", "electronic": "vivo.qmmm-electronic-forces"],
            outputs: [.init(name: "forces", kind: "vivo.qmmm-forces")],
            summary: "Identity-bound QM/MM interaction forces with MM reactions, fixed-distance link Jacobians and virtual-site redistribution.",
            configure: VivoPlatformOperations.empty, calculate: { _, inputs, budget in
                let document = try VivoPlatformOperations.input(VivoMolecularStructureDocument.self, "structure", inputs)
                let system = try VivoPlatformOperations.input(VivoClassicalSystem.self, "system", inputs)
                let frame = try VivoPlatformOperations.input(VivoMDQMMMPreparedFrame.self, "frame", inputs)
                let electronic = try VivoPlatformOperations.input(VivoQMMMElectronicForceResult.self, "electronic", inputs)
                let forces = try VivoQMMMForceMapper.assemble(document: document, system: system, frame: frame,
                    electronic: electronic, budget: budget)
                return ["forces": try VivoCanonicalJSON.encode(forces)]
            })]
    }
}
