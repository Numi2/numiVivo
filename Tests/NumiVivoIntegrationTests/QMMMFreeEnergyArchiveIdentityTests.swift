import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMFreeEnergyArchiveIdentityTests {
    private struct Fixture: Sendable {
        let system: VivoClassicalSystem
        let provider: VivoMDCandidateForceProvider
        let request: VivoQMMMFreeEnergyExecutionRequest
    }

    // Serialized archive shape, independent of the runner's private cursor type.
    private struct Cursor: Codable {
        let schema: String
        let numericalContract: String?
        let requestFingerprint: VivoFingerprint
        let completedWindows: Int
        let currentWindowProductionSteps: UInt64
        let currentMDCheckpoint: VivoFingerprint?
        let traces: [VivoQMMMUmbrellaTrace]
        let currentCoordinates: [Double]
        let currentEnergies: [Double]
    }

    private func fixture() throws -> Fixture {
        let structureID = try VivoCanonicalJSON.fingerprint(Data("archive-identity-host-fixture".utf8))
        let system = VivoClassicalSystem(identifier: "archive-identity-host-fixture", structureFingerprint: structureID,
            particles: [
                .init(index: 0, atomIndex: 0, typeIdentifier: "A", massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0),
                .init(index: 1, atomIndex: 1, typeIdentifier: "B", massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
            ])
        let systemID = try system.fingerprint()
        let providerID = try VivoCanonicalJSON.fingerprint(Data("unused-archive-identity-provider".utf8))
        let provider = try VivoMDCandidateForceProvider(fingerprint: providerID, retainedSystemFingerprint: systemID,
            boundary: .finiteCluster, supportsCellMoves: false, maximumAcceptedResidual: 1e-6,
            molecularConnectivitySystem: system) { _ in
                throw VivoChemistryError.invalid("host-only archive identity test must not evaluate a force provider")
            }
        let initial = VivoClassicalInitialState(systemFingerprint: systemID,
            positionsNM: [.zero, .init(0.2, 0, 0)])
        let windows = [0.2, 0.3].enumerated().map {
            VivoQMMMUmbrellaWindow(identifier: "window-\($0.offset)", centerNM: $0.element,
                forceConstantKJPerMolNM2: 100)
        }
        let sampling = VivoQMMMFreeEnergyRunRequest(
            coordinate: .init(identifier: "distance", kind: .distance, atomIndices: [0, 1]),
            windows: windows, initialStates: [initial, initial], randomSeeds: [17, 29],
            dynamics: .init(electrostatics: .cutoff), equilibrationSteps: 1,
            productionSteps: 2, sampleEvery: 1,
            analysis: .init(reactantRangeNM: 0.1 ... 0.2, dividingSurfaceNM: 0.3))
        return .init(system: system, provider: provider,
            request: .init(sampling: sampling, systemFingerprint: systemID, baseProviderFingerprint: providerID,
                execution: .init(productionStepsPerCheckpointBlock: 2)))
    }

    private func directory() throws -> URL {
        let root = FileManager.default.temporaryDirectory
            .appendingPathComponent("numivivo-qmmm-archive-identity-\(UUID().uuidString)")
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        return root
    }

    /// Pre-cancellation prevents native allocation even if an identity guard
    /// regresses. The expected identity error must precede the cancellation receipt.
    private func attempt(_ request: VivoQMMMFreeEnergyExecutionRequest, fixture: Fixture,
                         store: VivoArtifactStore, resume: VivoFingerprint? = nil) async throws
        -> VivoQMMMFreeEnergyExecutionReceipt {
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try await VivoQMMMFreeEnergyArchiveRunner.run(request, system: fixture.system,
                baseProvider: fixture.provider, store: store, resumeFrom: resume)
        }
        return try await task.value
    }

    @Test func legacyRequestsDecodeForInspectionButCannotExecute() async throws {
        let fixture = try fixture()
        // The integration target advances to v3; both earlier contracts must
        // remain inspectable while being excluded from execution.
        let contracts: [String?] = [nil, "numivivo.org/md-metal-numerics/v1", "numivivo.org/md-metal-numerics/v2"]
        try #require(contracts.allSatisfy { $0 != VivoMDExecutionIdentity.current })
        for contract in contracts {
            let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
            let store = try VivoArtifactStore(rootURL: root)
            var request = fixture.request
            request.numericalContract = contract
            let data = try VivoCanonicalJSON.encode(request)
            let decoded = try VivoCanonicalJSON.decode(VivoQMMMFreeEnergyExecutionRequest.self, from: data)
            #expect(decoded.numericalContract == contract)
            if contract == nil {
                let object = try #require(JSONSerialization.jsonObject(with: data) as? [String: Any])
                #expect(object["numericalContract"] == nil)
            }
            do {
                _ = try await attempt(decoded, fixture: fixture, store: store)
                Issue.record("old or absent request contract reached execution")
            } catch VivoArtifactValidationError.incompatible(let reason) {
                #expect(reason.contains("request numerical contract"))
            }
            #expect(try await store.list().isEmpty)
        }
    }

    @Test(arguments: [1, 2])
    func legacyWindowBoundariesCannotMixOrReanalyzeOldTraces(completedWindows: Int) async throws {
        let fixture = try fixture()
        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(fixture.request))
        let contracts: [String?] = [nil, "numivivo.org/md-metal-numerics/v1", "numivivo.org/md-metal-numerics/v2"]
        try #require(contracts.allSatisfy { $0 != VivoMDExecutionIdentity.current })
        for contract in contracts {
            let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
            let store = try VivoArtifactStore(rootURL: root)
            // These are declared synthetic identity fixtures, not sampled PMFs.
            let traces = fixture.request.sampling.windows.prefix(completedWindows).enumerated().map {
                VivoQMMMUmbrellaTrace(window: $0.element, randomSeed: fixture.request.sampling.randomSeeds[$0.offset],
                    coordinateNM: [$0.element.centerNM, $0.element.centerNM + 0.001],
                    potentialEnergyKJPerMol: [0, 0.1])
            }
            // A current request hash isolates the cursor contract guard. Both a
            // completed-window boundary and an all-complete cursor lack MD state.
            let cursor = Cursor(schema: "numivivo.org/qmmm-free-energy-execution-checkpoint/v1",
                numericalContract: contract, requestFingerprint: requestID, completedWindows: completedWindows,
                currentWindowProductionSteps: 0, currentMDCheckpoint: nil, traces: traces,
                currentCoordinates: [], currentEnergies: [])
            let data = try VivoCanonicalJSON.encode(cursor)
            #expect(try VivoCanonicalJSON.decode(Cursor.self, from: data).numericalContract == contract)
            let stored = try await store.put(data: data, kind: "qmmm-free-energy-execution-checkpoint", mediaType: "application/json")
            do {
                _ = try await attempt(fixture.request, fixture: fixture, store: store, resume: stored.fingerprint)
                Issue.record("old window prefix reached continuation or analysis")
            } catch VivoArtifactValidationError.incompatible(let reason) {
                #expect(reason.contains("cursor numerical contract"))
            }
            #expect(try await store.data(for: stored.fingerprint) == data)
        }
    }

    @Test func freshExecutionPublishesCurrentContractBeforeCancellation() async throws {
        let fixture = try fixture()
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let store = try VivoArtifactStore(rootURL: root)
        #expect(fixture.request.numericalContract == VivoMDExecutionIdentity.current)
        let receipt = try await attempt(fixture.request, fixture: fixture, store: store)
        #expect(receipt.status == .cancelled)
        let data = try await store.data(for: receipt.checkpoint)
        let cursor = try VivoCanonicalJSON.decode(Cursor.self, from: data)
        #expect(cursor.numericalContract == VivoMDExecutionIdentity.current)
        #expect(cursor.requestFingerprint == receipt.requestFingerprint)
        #expect(cursor.completedWindows == 0 && cursor.currentMDCheckpoint == nil && cursor.traces.isEmpty)
        let resumed = try await attempt(fixture.request, fixture: fixture, store: store, resume: receipt.checkpoint)
        #expect(resumed.status == .cancelled && resumed.checkpoint == receipt.checkpoint)
    }
}
