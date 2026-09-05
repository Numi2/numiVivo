import Foundation
import NumiVivoKit

private struct CheckFailure: Error { let message: String }

private actor SuspendedBackend: VivoSurrogateBackend {
    nonisolated let backendID = "numivivo.checks.suspended"
    nonisolated let modelFingerprint = "synthetic-model"
    nonisolated let inputWidth = 1
    nonisolated let outputWidth = 1
    private var continuation: CheckedContinuation<Void, Never>?
    private var observers: [CheckedContinuation<Void, Never>] = []
    func predict(_ batch: VivoSurrogateBatch) async throws -> VivoSurrogatePrediction {
        await withCheckedContinuation { handle in
            continuation = handle
            for observer in observers { observer.resume() }
            observers = []
        }
        return try .init(batchSize: batch.batchSize, outputWidth: 1, values: batch.values.map { 2 * $0 },
                         normalizedUncertainty: [Float](repeating: 0.01, count: batch.batchSize),
                         backend: backendID, modelFingerprint: modelFingerprint)
    }
    func waitUntilBlocked() async {
        if continuation != nil { return }
        await withCheckedContinuation { observers.append($0) }
    }
    func release() { continuation?.resume(); continuation = nil }
}

/// Apple integration checks, not a claim that they were run during source authoring.
@main struct TargetEngagementChecks {
    static func main() async throws {
        var checks = 0
        func check(_ condition: Bool, _ name: String) throws {
            guard condition else { throw CheckFailure(message: name) }
            checks += 1
        }
        let example = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
            .deletingLastPathComponent().deletingLastPathComponent()
        let experiment = try VivoKineticsDocumentIO.read(VivoTargetEngagementExperiment.self,
            from: example.appendingPathComponent("synthetic-pulse.json"))
        let reference = try VivoTargetEngagementRunRecord.reference(experiment)
        try reference.validate(); checks += 1
        let compiled = try VivoTargetEngagementCompiler.compile(experiment)
        try check(compiled.targetSpeciesIndices.count == 4, "compiled target index mapping")
        try VivoProgramExecutionContract.validate(pack: compiled.pack,
            configuration: .init(fidelity: compiled.pack.header.fidelity, maximumSubsteps: 1))
        checks += 1
        do {
            try VivoProgramExecutionContract.validate(pack: compiled.pack,
                configuration: .init(fidelity: .stochastic, maximumSubsteps: 2))
            throw CheckFailure(message: "legacy stochastic resampling accepted")
        } catch VivoRuntimeError.invalidConfiguration(let message) {
            try check(message.contains("maximumSubsteps"), "legacy retry guard diagnostic")
        }
        do {
            try VivoCalibrationStrategy().validate(parameterCount: Int.max)
            throw CheckFailure(message: "calibration count overflow accepted")
        } catch VivoArtifactValidationError.invalid { checks += 1 }
        let feature = VivoSurrogateFeature(id: "x", unit: "fraction", minimum: 0, maximum: 10)
        func contract(interval: UInt32 = 1, uncertainty: VivoSurrogateUncertaintyKind = .predictedStandardDeviation) throws -> VivoSurrogateContract {
            try .init(id: "synthetic-contract", modelFingerprint: "synthetic-model",
                trainingDataFingerprint: "synthetic-data", mechanismPackFingerprint: "synthetic-mechanism",
                inputs: [feature], outputs: [feature], uncertaintyKind: uncertainty,
                maximumNormalizedUncertainty: uncertainty == .none ? 0 : 0.1,
                maximumConsecutiveAcceptedSteps: interval, mandatoryAuthorityInterval: interval)
        }
        let affine = try VivoAffineSurrogateBackend(modelFingerprint: "synthetic-model", inputWidth: 1, outputWidth: 1,
            weights: [2], bias: [0], outputUncertainty: [0.01])
        let gate = try VivoSurrogateAuthorityGate(contract: contract(), backend: affine)
        let batch = try VivoSurrogateBatch(batchSize: 1, inputWidth: 1, values: [1])
        let first = await gate.evaluate(batch)
        try check(first.requiresAuthoritativeEvaluation && !first.accepted, "initial refresh required")
        guard let request = first.authorityRequest else { throw CheckFailure(message: "missing refresh identity") }
        let repeated = await gate.evaluate(batch)
        try check(repeated.authorityRequest == request, "refresh remains pending")
        let resultID = try VivoCanonicalJSON.fingerprint(Data("synthetic authoritative result: 2".utf8)).hex
        do {
            try await gate.recordAuthoritativeEvaluation(request: request, resultFingerprint: resultID,
                values: [2], numericalChecksPassed: false)
            throw CheckFailure(message: "failed authoritative result accepted")
        } catch VivoSurrogateError.invalidAuthorityReceipt { checks += 1 }
        try check(await gate.pendingAuthorityRequest() == request, "failed receipt preserves request")
        try await gate.recordAuthoritativeEvaluation(request: request, resultFingerprint: resultID,
            values: [2], numericalChecksPassed: true)
        let accepted = await gate.evaluate(batch)
        try check(accepted.accepted, "qualified in-domain surrogate step")
        let refresh = await gate.evaluate(batch)
        try check(!refresh.accepted && refresh.authorityRequest?.generation == request.generation + 1,
                  "refresh generation advances only on authoritative receipt")
        do {
            try await gate.recordAuthoritativeEvaluation(request: request, resultFingerprint: resultID,
                values: [2], numericalChecksPassed: true)
            throw CheckFailure(message: "stale receipt accepted")
        } catch VivoSurrogateError.invalidAuthorityReceipt { checks += 1 }
        let noUncertainty = try VivoSurrogateAuthorityGate(contract: contract(uncertainty: .none), backend: affine)
        let noFirst = await noUncertainty.evaluate(batch)
        try await noUncertainty.recordAuthoritativeEvaluation(request: noFirst.authorityRequest!, resultFingerprint: resultID,
            values: [2], numericalChecksPassed: true)
        let rejectedUnknown = await noUncertainty.evaluate(batch)
        try check(rejectedUnknown.reason == .unavailableUncertainty, "unknown uncertainty never accepts")

        let suspended = SuspendedBackend()
        let serialGate = try VivoSurrogateAuthorityGate(contract: contract(interval: 8), backend: suspended)
        let initialize = await serialGate.evaluate(batch)
        try await serialGate.recordAuthoritativeEvaluation(request: initialize.authorityRequest!, resultFingerprint: resultID,
            values: [2], numericalChecksPassed: true)
        async let pending = serialGate.evaluate(batch)
        await suspended.waitUntilBlocked()
        let busy = await serialGate.evaluate(batch)
        try check(busy.reason == .evaluationInFlight && !busy.accepted, "gate actor single-flight across await")
        await suspended.release()
        let completed = await pending
        try check(completed.accepted, "reserved prediction completes")

        #if canImport(CoreML)
        do {
            _ = try VivoCoreMLSurrogateBackend(modelURL: URL(fileURLWithPath: "/nonexistent/synthetic.mlmodelc"),
                modelFingerprint: "synthetic-model", inputWidth: 1, outputWidth: 1)
            throw CheckFailure(message: "missing Core ML uncertainty accepted")
        } catch VivoSurrogateError.backendContractMismatch { checks += 1 }
        #endif

        let temporary = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString, isDirectory: true)
        defer { try? FileManager.default.removeItem(at: temporary) }
        let document = temporary.appendingPathComponent("experiment.json")
        try VivoKineticsDocumentIO.write(VivoCanonicalJSON.encode(experiment), to: document)
        let restored = try VivoKineticsDocumentIO.read(VivoTargetEngagementExperiment.self, from: document)
        try check(restored == experiment, "bounded document roundtrip")
        do {
            try VivoKineticsDocumentIO.write(Data("replacement".utf8), to: document)
            throw CheckFailure(message: "existing output overwritten without permission")
        } catch VivoKineticsError.invalid { checks += 1 }

        if CommandLine.arguments.contains("--metal") {
            let gpu = try await VivoTargetEngagementMetalRunner.run([experiment])
            try check(gpu.results[0].samples.count == reference.result.samples.count, "Metal output shape")
            for (actual, expected) in zip(gpu.results[0].samples, reference.result.samples) {
                try check(actual.timeSeconds == expected.timeSeconds, "Metal exact observation clock")
                try check(abs(actual.drugOccupancy - expected.drugOccupancy) < 1e-3, "Metal/reference occupancy agreement")
            }
        }
        print("\(checks) Apple integration checks passed; Metal executed: \(CommandLine.arguments.contains("--metal")).")
    }
}
