import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct PlatformQMMMWorkflowTests {
    @Test func registeredCheckpointAndForceOperationsReconstructTheirOutputs() async throws {
        let carbon = try #require(VivoElement.from(symbol: "C"))
        let atoms = (0..<3).map { VivoMolecularAtom(index: UInt32($0), name: "C\($0)", element: carbon) }
        let cell = VivoPeriodicCell(a: .init(2,0,0), b: .init(0,2,0), c: .init(0,0,2))
        let atomPositions: [VivoVector3D] = [.init(0.125,0,0),.init(1.875,0,0),.init(0.75,0.5,0.25)]
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "workflow-qmmm-fixture", atoms: atoms,
            bonds: [.init(atomA: 0, atomB: 1)], conformers: [.init(positionsNM: atomPositions)], periodicCell: cell))
        let order = [2,0,1]
        let particles = order.enumerated().map { slot, atom in VivoClassicalParticle(index: UInt32(slot),
            atomIndex: UInt32(atom), typeIdentifier: "C", massDa: 12, chargeE: atom == 2 ? 0.25 : 0,
            sigmaNM: 0, epsilonKJPerMol: 0) }
        let system = VivoClassicalSystem(identifier: "workflow-qmmm-fixture", structureFingerprint: document.structureFingerprint,
            particles: particles)
        let identity = try system.fingerprint()
        let checkpoint = VivoMDCheckpoint(systemFingerprint: identity, configurationFingerprint: identity,
            acceptedStep: 17, timePS: 0.034, positionsNM: order.map { atomPositions[$0] },
            velocitiesNMPerPS: [.zero,.zero,.zero], periodicCell: cell)
        let request = VivoQMMMRegionRequest(qmAtomIndices: [0], alphaElectrons: 4, betaElectrons: 3)
        let configuration = try VivoPlatformOperations.json(VivoWorkflowMDQMMMConfiguration(request: request))
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: identity)
        for identifier in ["vivo.platform.qmmm-free-energy-analyze","vivo.platform.qmmm-free-energy-rate",
                           "vivo.platform.qmmm-transmission-analyze","vivo.platform.qmmm-apply-transmission",
                           "vivo.platform.qmmm-replicated-free-energy-rate","vivo.platform.qmmm-apply-replicated-rate",
                           "vivo.platform.qmmm-chemical-qualification","vivo.platform.qmmm-chemical-state-network",
                           "vivo.platform.qmmm-chemical-exchange-network"] {
            _ = try registry.definition(identifier)
        }
        let prepare = try registry.definition("vivo.platform.md-qmmm").operation
        let inputs = ["structure": try document.canonicalData(), "system": try VivoCanonicalJSON.encode(system),
                      "checkpoint": try VivoCanonicalJSON.encode(checkpoint)]
        let outputs = try await prepare.execute(configuration, inputs, .init())
        try prepare.validateOutputs(configuration, inputs, outputs, .init())
        let frame = try VivoPlatformOperations.input(VivoMDQMMMPreparedFrame.self, "frame", outputs)
        let electronic = try VivoPlatformOperations.input(VivoElectronicSystem.self, "system", outputs)
        let mapping = try VivoPlatformOperations.input(VivoQMMMFiniteCluster.self, "mapping", outputs)
        #expect(mapping == frame.cluster && mapping.atomToParticle == [1,2,0])
        #expect(electronic == frame.region.electronicSystem && frame.sourceStep == 17)
        #expect(frame.sourceCheckpointFingerprint == (try checkpoint.fingerprint()))
        #expect(frame.region.links.count == 1 && electronic.pointCharges.map(\.classicalParticleIndex) == [0])
        #expect(electronic.nuclei.map(\.structureAtomIndex) == [0,nil])
        let direct = try VivoQMMMCompiler.prepareMD(document: document, system: system, checkpoint: checkpoint, request: request)
        #expect(frame == direct)

        // A center-distance fixture exercises the registered force adapter without
        // asserting chemical validation or depending on a particular paper/backend.
        let centers = electronic.nuclei.map(\.positionBohr) + electronic.pointCharges.map(\.positionBohr)
        var energy = 0.0, forces = [SIMD3<Double>](repeating: .zero, count: centers.count)
        for i in centers.indices { for j in 0..<i {
            let delta = centers[i] - centers[j]
            energy += 0.0005 * (delta.x*delta.x + delta.y*delta.y + delta.z*delta.z)
            forces[i] -= delta * 0.001; forces[j] += delta * 0.001
        } }
        let vectors = forces.map { VivoVector3D($0.x,$0.y,$0.z) }
        let evaluated = try VivoQMMMElectronicForceResult(system: electronic, energyHartree: energy,
            nucleusForcesHartreePerBohr: Array(vectors.prefix(electronic.nuclei.count)),
            pointChargeForcesHartreePerBohr: Array(vectors.dropFirst(electronic.nuclei.count)),
            derivativeMethod: "synthetic-workflow-center-fixture")
        let assemble = try registry.definition("vivo.platform.qmmm-forces").operation
        let forceInputs = ["structure": inputs["structure"]!, "system": inputs["system"]!,
                           "frame": outputs["frame"]!, "electronic": try VivoCanonicalJSON.encode(evaluated)]
        let result = try await assemble.execute(.object([:]), forceInputs, .init())
        try assemble.validateOutputs(.object([:]), forceInputs, result, .init())
        let mapped = try VivoPlatformOperations.input(VivoQMMMForceResult.self, "forces", result)
        #expect(mapped.atomToParticle == [1,2,0] && mapped.particleForcesHartreePerBohr.count == 3)
        #expect(mapped.sourceCheckpointFingerprint == frame.sourceCheckpointFingerprint)
        #expect(mapped.netForceHartreePerBohr.norm < 1e-12 && mapped.netTorqueHartree.norm < 1e-12)
        var altered = outputs; altered["mapping"] = Data("null".utf8)
        #expect(throws: (any Error).self) { try prepare.validateOutputs(configuration, inputs, altered, .init()) }
    }
}
