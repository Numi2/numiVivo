import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct PlatformSnapshotTests {
    private func identity() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data("platform-snapshot-tests".utf8))
    }
    private func directory() throws -> URL {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-snapshot-\(UUID().uuidString)")
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        return root
    }
    private func record<T: Encodable>(_ value: T, _ name: String) throws {
        guard let path = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let root = URL(fileURLWithPath: path)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        try VivoCanonicalJSON.encode(value).write(to: root.appendingPathComponent(name + ".json"), options: .atomic)
    }
    @Test func acceptedMDSnapshotFeedsElectronicWorkflowWithoutFileConversion() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: identity())
        let store = try VivoArtifactStore(rootURL: root), recipe = try VivoWorkflowTemplates.mdElectronicAnalysis()
        let executor = VivoWorkflowExecutor(store: store, registry: registry), workflow = VivoChemistryWorkflow(store: store)
        let result = try await executor.run(recipe)
        try record(result.report, "platform-md-electronic-analysis")
        try #require(result.report.allTasksSucceeded)
        #expect(result.report.nodes.count == 9)
        let mapArtifact = try #require(result.report.exports.first(where: { $0.name == "snapshot-mapping" })?.artifact)
        let mapping = try VivoCanonicalJSON.decode(VivoWorkflowSnapshotMapping.self,
            from: await workflow.payload(artifact: mapArtifact.artifact, expectedKind: mapArtifact.kind))
        #expect(mapping.atomToParticle == [0, 1] && mapping.step == 20)
        #expect(mapping.snapshotStructure != mapping.sourceStructure)
        let cpArtifact = try #require(result.report.exports.first(where: { $0.name == "checkpoint" })?.artifact)
        let checkpoint = try VivoCanonicalJSON.decode(VivoMDCheckpoint.self,
            from: await workflow.payload(artifact: cpArtifact.artifact, expectedKind: cpArtifact.kind))
        let systemArtifact = try #require(result.report.nodes.first(where: { $0.identifier == "system" })?.outcome.outputs?.first)
        let electronic = try VivoCanonicalJSON.decode(VivoElectronicSystem.self,
            from: await workflow.payload(artifact: systemArtifact.artifact, expectedKind: systemArtifact.kind))
        #expect(electronic.nuclei.map(\.structureAtomIndex) == [0, 1])
        for (i, atom) in electronic.nuclei.enumerated() {
            let p = checkpoint.positionsNM[Int(mapping.atomToParticle[i])]
            #expect(atom.positionBohr == SIMD3<Double>(p.x,p.y,p.z) / VivoAtomicUnits.bohrInNM)
        }
        let cached = try await executor.run(recipe)
        #expect(cached.report.nodes.allSatisfy { $0.outcome.reused })
    }

    @Test func snapshotMappingRetainsPeriodicityAndDoesNotAssumeParticleOrder() async throws {
        let recipe = try VivoWorkflowTemplates.molecularDynamics()
        func inline<T: Decodable>(_ name: String, _ type: T.Type) throws -> T {
            guard let artifact = recipe.artifacts.first(where: { $0.identifier == name }),
                  case .json(_, let payload) = artifact.source else { throw VivoChemistryError.invalid("missing test input") }
            return try VivoPlatformOperations.decode(type, payload)
        }
        let document = try inline("structure", VivoMolecularStructureDocument.self)
        var system = try inline("system", VivoClassicalSystem.self)
        system.particles[0].atomIndex = 1; system.particles[1].atomIndex = 0
        let cell = VivoPeriodicCell(a: .init(4,0,0), b: .init(0,4,0), c: .init(0,0,4))
        var checkpoint = VivoMDCheckpoint(systemFingerprint: try system.fingerprint(), configurationFingerprint: try identity(),
            acceptedStep: 7, timePS: 0.5, positionsNM: [.init(1,0,0), .init(2,0,0)],
            velocitiesNMPerPS: [.zero, .zero], periodicCell: cell)
        let operation = try VivoPlatformSnapshotOperations.definition(implementationFingerprint: identity()).operation
        var inputs = ["structure": try document.canonicalData(), "system": try VivoCanonicalJSON.encode(system),
                      "checkpoint": try VivoCanonicalJSON.encode(checkpoint)]
        let outputs = try await operation.execute(.object([:]), inputs, .init())
        try operation.validateOutputs(.object([:]), inputs, outputs, .init())
        let snapshot = try VivoPlatformOperations.input(VivoMolecularStructureDocument.self, "structure", outputs)
        let mapping = try VivoPlatformOperations.input(VivoWorkflowSnapshotMapping.self, "mapping", outputs)
        #expect(mapping.atomToParticle == [1, 0])
        #expect(snapshot.structure.conformers[0].positionsNM == [.init(2,0,0), .init(1,0,0)])
        #expect(snapshot.structure.periodicCell == cell)
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: identity())
        let electronic = try registry.definition("vivo.platform.structure-electronic-system").operation
        do {
            _ = try await electronic.execute(VivoPlatformOperations.json(VivoWorkflowElectronicSystem(alphaElectrons: 1, betaElectrons: 1)),
                ["structure": outputs["structure"]!], .init())
            Issue.record("periodic snapshot was silently treated as isolated")
        } catch {}
        system.particles[1].atomIndex = 1
        checkpoint.systemFingerprint = try system.fingerprint()
        inputs["system"] = try VivoCanonicalJSON.encode(system)
        inputs["checkpoint"] = try VivoCanonicalJSON.encode(checkpoint)
        do { _ = try await operation.execute(.object([:]), inputs, .init()); Issue.record("duplicate snapshot ownership accepted") } catch {}
    }
    @Test func structureImportNeverDropsAdditionalMoleculeRecords() async throws {
        let recipe = try VivoWorkflowTemplates.molecularAnalysis()
        guard case .json(_, let payload) = recipe.artifacts[0].source else { throw VivoChemistryError.invalid("test source") }
        let document = try VivoPlatformOperations.decode(VivoMolecularStructureDocument.self, payload)
        let sdf = try VivoSDF.write(document.structure)
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: identity())
        let operation = try registry.definition("vivo.platform.structure-import").operation
        let cfg = try VivoPlatformOperations.json(VivoWorkflowStructureImport(format: .sdf))
        let single = ["source": try VivoCanonicalJSON.encode(sdf)]
        let outputs = try await operation.execute(cfg, single, .init())
        try operation.validateOutputs(cfg, single, outputs, .init())
        let result = try VivoPlatformOperations.input(VivoMolecularStructureDocument.self, "structure", outputs)
        #expect(result.structure.atoms.count == 2)
        // A missing final separator must not let the second molecule disappear.
        for text in [sdf + sdf, sdf + sdf.replacingOccurrences(of: "$$$$", with: "")] {
            do {
                _ = try await operation.execute(cfg, ["source": VivoCanonicalJSON.encode(text)], .init())
                Issue.record("additional SDF record silently dropped")
            } catch {}
        }
        for (format, text) in [(VivoStructureFormat.smiles, "C first\nO second\n"),
                               (.mol2, "@<TRIPOS>MOLECULE\nfirst\n@<TRIPOS>MOLECULE\nsecond\n")] {
            do {
                _ = try await operation.execute(VivoPlatformOperations.json(VivoWorkflowStructureImport(format: format)),
                    ["source": VivoCanonicalJSON.encode(text)], .init())
                Issue.record("additional molecule record silently dropped")
            } catch {}
        }
    }
}
