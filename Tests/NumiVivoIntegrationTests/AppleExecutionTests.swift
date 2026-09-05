import Foundation
@preconcurrency import Metal
import Testing
@testable import NumiVivoKit
import NumiVivoShaders

/// A real unified-memory GPU is required. Absent hardware fails qualification.
@Suite(.serialized) struct AppleExecutionTests {
    private func record<T: Encodable>(_ value: T, _ name: String) throws {
        guard let directory = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let url = URL(fileURLWithPath: directory, isDirectory: true)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: true)
        try VivoCanonicalJSON.encode(value).write(to: url.appendingPathComponent(name + ".json"), options: .atomic)
    }
    private func device() throws -> MTLDevice {
        let device = try VivoMetalDeviceSelector.productionDevice()
        #expect(device.hasUnifiedMemory)
        return device
    }
    private func experiment() throws -> VivoTargetEngagementExperiment {
        let root = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
            .deletingLastPathComponent().deletingLastPathComponent()
        return try VivoCanonicalJSON.decode(VivoTargetEngagementExperiment.self,
            from: Data(contentsOf: root.appendingPathComponent("Examples/target-engagement/synthetic-pulse.json")))
    }
    private func rejects(_ operation: () async throws -> Void) async {
        do { try await operation(); Issue.record("Expected an explicit operation rejection") }
        catch { }
    }
    @Test func completePipelineCompilation() async throws {
        let device = try device(), catalog = try NumiVivoPipelineCatalog(device: device)
        try await catalog.preloadAll()
        struct Report: Encodable { let device: String; let registryID: UInt64; let kernels: [String] }
        try record(Report(device: device.name, registryID: device.registryID,
            kernels: NumiVivoKernel.allCases.map(\.rawValue)), "pipeline-compilation")
    }
    @Test func molecularTransactionsAndInterventions() async throws {
        let compiled = try VivoTargetEngagementCompiler.compile(experiment())
        let configuration = VivoRuntimeConfiguration(fidelity: .init(rawValue: 1)!, environmentCount: 2,
            timeStep: 0.001, minimumTimeStep: 0.0001, maximumTimeStep: 0.01)
        let runtime = try await VivoTransactionalMolecularRuntime.make(pack: compiled.pack,
            configuration: configuration, device: device())
        let initial = try await runtime.resumeCheckpoint()
        try await runtime.apply(intervention: .coupling([
            .init(speciesIndex: compiled.drugSpeciesIndex, laneIndex: 0, value: 1e-6),
            .init(speciesIndex: compiled.drugSpeciesIndex, laneIndex: 1, value: 2e-6)]))
        let changed = try await runtime.snapshot()
        #expect(changed.absoluteTime == 0 && changed.stepIndex == 0)
        #expect(changed.value(species: compiled.drugSpeciesIndex, lane: 0) == Float(1e-6))
        #expect(changed.value(species: compiled.drugSpeciesIndex, lane: 1) == Float(2e-6))
        await rejects { try await runtime.apply(intervention: .coupling([
            .init(speciesIndex: compiled.targetSpeciesIndices[0], laneIndex: 0, value: 0)])) }
        #expect(try await runtime.snapshot() == changed)
        let pending = try await runtime.prepareStep(.init(timeStep: 0.001, permitAdaptiveReduction: false))
        try #require(pending.canCommit)
        await rejects { _ = try await runtime.snapshot() }
        await rejects { try await runtime.apply(intervention: .coupling([
            .init(speciesIndex: compiled.drugSpeciesIndex, laneIndex: 0, value: 0)])) }
        try await runtime.discardPreparedStep(transactionID: pending.transactionID)
        #expect(try await runtime.snapshot() == changed)
        let committed = try await runtime.step(.init(timeStep: 0.001, permitAdaptiveReduction: false))
        #expect(committed.certificate.committed)
        let accepted = try await runtime.resumeCheckpoint()
        try await runtime.restore(initial)
        #expect(try await runtime.snapshot().stepIndex == 0)
        try await runtime.restore(accepted)
        #expect(try await runtime.snapshot().stepIndex == 1)
        try await runtime.apply(intervention: .reversibleShutdown(reason: "integration-check"))
        await rejects { _ = try await runtime.step() }
        try await runtime.resume()
        #expect(try await runtime.step().certificate.committed)
        try await runtime.apply(intervention: .permanentShutdown(reason: "integration-check"))
        await rejects { try await runtime.resume() }
        await rejects { try await runtime.restore(initial) }
        try record(accepted, "molecular-accepted-checkpoint")
    }
    @Test func targetKineticsOnMetalAgainstFP64() async throws {
        let input = try experiment()
        let reference = try VivoTargetEngagementReference.run(input)
        let result = try await VivoTargetEngagementMetalRunner.run([input],
            policy: .init(maximumTimeStepSeconds: 1.0 / 128.0), device: device())
        #expect(result.results.count == 1)
        let native = try #require(result.results.first)
        #expect(native.samples.count == reference.samples.count)
        var maximum = 0.0
        for (a,b) in zip(native.samples, reference.samples) {
            maximum = max(maximum, abs(a.drugOccupancy-b.drugOccupancy), abs(a.covalentOccupancy-b.covalentOccupancy))
        }
        #expect(maximum < 2e-4)
        try record(result, "target-metal-result")
        try record(reference, "target-fp64-reference")
    }
    @Test func harmonicMDOnMetalAndExactRestart() async throws {
        let identity = try VivoCanonicalJSON.fingerprint(Data("explicit two-atom harmonic regression".utf8))
        let particles = (0..<2).map { VivoClassicalParticle(index: UInt32($0), atomIndex: UInt32($0),
            typeIdentifier: "regression-H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0) }
        let system = VivoClassicalSystem(identifier: "harmonic-pair", structureFingerprint: identity,
            particles: particles, bonds: [.init(a: 0, b: 1, lengthNM: 0.1, forceConstant: 1000)])
        let configuration = VivoMDConfiguration(timeStepPS: 0.0001, electrostatics: .cutoff,
            ensemble: .nve, thermostat: .none, targetTemperatureK: nil, frictionPerPS: nil,
            neighborListEnabled: false)
        let initial = VivoClassicalInitialState(systemFingerprint: try system.fingerprint(),
            positionsNM: [.init(-0.06,0,0),.init(0.06,0,0)], periodicCell: nil)
        let runtime = try await VivoMDMetalRuntime.make(system: system, initialState: initial,
            configuration: configuration, device: device())
        let before = try await runtime.observables()
        #expect(abs(before.potentialEnergyKJPerMol-0.2) < 2e-6)
        for _ in 0..<100 { #expect(try await runtime.step().committed) }
        let after = try await runtime.observables(), snapshot = try await runtime.snapshot()
        #expect(abs(after.totalEnergyKJPerMol-before.totalEnergyKJPerMol) < 2e-5)
        #expect((snapshot.positionsNM[0]+snapshot.positionsNM[1]).norm < 1e-7)
        let checkpoint = try await runtime.checkpoint()
        let resumed = try await VivoMDMetalRuntime.restore(system: system, configuration: configuration,
            checkpoint: checkpoint, device: device())
        for _ in 0..<10 {
            #expect(try await runtime.step().committed)
            #expect(try await resumed.step().committed)
        }
        let uninterrupted = try await runtime.checkpoint(), restarted = try await resumed.checkpoint()
        #expect(uninterrupted == restarted)
        try record(before, "md-harmonic-initial")
        try record(after, "md-harmonic-after100")
        try record(checkpoint, "md-harmonic-checkpoint")
    }
}
