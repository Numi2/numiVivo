import CryptoKit
import Foundation

public struct VivoSurrogateFeature: Codable, Hashable, Sendable {
    public let id: String
    public let unit: String
    public let minimum: Float
    public let maximum: Float
    public let normalizationOffset: Float
    public let normalizationScale: Float

    public init(id: String, unit: String, minimum: Float, maximum: Float,
                normalizationOffset: Float = 0, normalizationScale: Float = 1) {
        self.id = id; self.unit = unit; self.minimum = minimum; self.maximum = maximum
        self.normalizationOffset = normalizationOffset; self.normalizationScale = normalizationScale
    }
    public func normalize(_ value: Float) -> Float { (value - normalizationOffset) * normalizationScale }
    public func denormalize(_ value: Float) -> Float { value / normalizationScale + normalizationOffset }
}

public enum VivoSurrogateUncertaintyKind: String, Codable, Sendable {
    case none, predictedStandardDeviation, ensembleVariance, conformalRadius
}

/// A declared applicability domain, not a biological validation certificate.
public struct VivoSurrogateContract: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let id: String
    public let modelFingerprint: String
    public let trainingDataFingerprint: String
    public let mechanismPackFingerprint: String
    public let hostContextFingerprint: String?
    public let inputs: [VivoSurrogateFeature]
    public let outputs: [VivoSurrogateFeature]
    public let uncertaintyKind: VivoSurrogateUncertaintyKind
    public let maximumNormalizedUncertainty: Float
    public let maximumNormalizedExtrapolation: Float
    public let maximumConsecutiveAcceptedSteps: UInt32
    public let mandatoryAuthorityInterval: UInt32
    public let fingerprint: String

    public init(id: String, modelFingerprint: String, trainingDataFingerprint: String,
                mechanismPackFingerprint: String, hostContextFingerprint: String? = nil,
                inputs: [VivoSurrogateFeature], outputs: [VivoSurrogateFeature],
                uncertaintyKind: VivoSurrogateUncertaintyKind, maximumNormalizedUncertainty: Float,
                maximumNormalizedExtrapolation: Float = 0,
                maximumConsecutiveAcceptedSteps: UInt32 = 32, mandatoryAuthorityInterval: UInt32 = 32) throws {
        let unsigned = VivoSurrogateContract(schemaVersion: 1, id: id, modelFingerprint: modelFingerprint,
            trainingDataFingerprint: trainingDataFingerprint, mechanismPackFingerprint: mechanismPackFingerprint,
            hostContextFingerprint: hostContextFingerprint, inputs: inputs, outputs: outputs,
            uncertaintyKind: uncertaintyKind, maximumNormalizedUncertainty: maximumNormalizedUncertainty,
            maximumNormalizedExtrapolation: maximumNormalizedExtrapolation,
            maximumConsecutiveAcceptedSteps: maximumConsecutiveAcceptedSteps,
            mandatoryAuthorityInterval: mandatoryAuthorityInterval, fingerprint: "")
        try unsigned.validate()
        self = unsigned.withFingerprint(try Self.fingerprint(unsigned))
    }

    private init(schemaVersion: UInt32, id: String, modelFingerprint: String, trainingDataFingerprint: String,
                 mechanismPackFingerprint: String, hostContextFingerprint: String?,
                 inputs: [VivoSurrogateFeature], outputs: [VivoSurrogateFeature],
                 uncertaintyKind: VivoSurrogateUncertaintyKind, maximumNormalizedUncertainty: Float,
                 maximumNormalizedExtrapolation: Float, maximumConsecutiveAcceptedSteps: UInt32,
                 mandatoryAuthorityInterval: UInt32, fingerprint: String) {
        self.schemaVersion = schemaVersion; self.id = id; self.modelFingerprint = modelFingerprint
        self.trainingDataFingerprint = trainingDataFingerprint; self.mechanismPackFingerprint = mechanismPackFingerprint
        self.hostContextFingerprint = hostContextFingerprint; self.inputs = inputs; self.outputs = outputs
        self.uncertaintyKind = uncertaintyKind; self.maximumNormalizedUncertainty = maximumNormalizedUncertainty
        self.maximumNormalizedExtrapolation = maximumNormalizedExtrapolation
        self.maximumConsecutiveAcceptedSteps = maximumConsecutiveAcceptedSteps
        self.mandatoryAuthorityInterval = mandatoryAuthorityInterval; self.fingerprint = fingerprint
    }

    private func withFingerprint(_ value: String) -> Self {
        .init(schemaVersion: schemaVersion, id: id, modelFingerprint: modelFingerprint,
              trainingDataFingerprint: trainingDataFingerprint, mechanismPackFingerprint: mechanismPackFingerprint,
              hostContextFingerprint: hostContextFingerprint, inputs: inputs, outputs: outputs,
              uncertaintyKind: uncertaintyKind, maximumNormalizedUncertainty: maximumNormalizedUncertainty,
              maximumNormalizedExtrapolation: maximumNormalizedExtrapolation,
              maximumConsecutiveAcceptedSteps: maximumConsecutiveAcceptedSteps,
              mandatoryAuthorityInterval: mandatoryAuthorityInterval, fingerprint: value)
    }

    public func validate() throws {
        guard schemaVersion == 1, !id.isEmpty, !modelFingerprint.isEmpty, !trainingDataFingerprint.isEmpty,
              !mechanismPackFingerprint.isEmpty, !inputs.isEmpty, !outputs.isEmpty,
              Set(inputs.map(\.id)).count == inputs.count, Set(outputs.map(\.id)).count == outputs.count,
              maximumNormalizedUncertainty.isFinite, maximumNormalizedUncertainty >= 0,
              maximumNormalizedExtrapolation.isFinite, maximumNormalizedExtrapolation >= 0,
              maximumConsecutiveAcceptedSteps > 0, mandatoryAuthorityInterval > 0,
              mandatoryAuthorityInterval <= maximumConsecutiveAcceptedSteps else {
            throw VivoSurrogateError.invalidContract("top-level contract is inconsistent")
        }
        for feature in inputs + outputs {
            guard !feature.id.isEmpty, !feature.unit.isEmpty, feature.minimum.isFinite, feature.maximum.isFinite,
                  feature.minimum <= feature.maximum, feature.normalizationOffset.isFinite,
                  feature.normalizationScale.isFinite, feature.normalizationScale != 0 else {
                throw VivoSurrogateError.invalidContract("feature \(feature.id) has invalid bounds or normalization")
            }
        }
        if uncertaintyKind == .none, maximumNormalizedUncertainty != 0 {
            throw VivoSurrogateError.invalidContract("uncertainty threshold must be zero when uncertaintyKind is none")
        }
    }

    /// A decoded contract must bind its actual fields, not only repeat a claimed hash.
    public func validateIntegrity() throws {
        try validate()
        guard fingerprint == (try Self.fingerprint(withFingerprint(""))) else {
            throw VivoSurrogateError.invalidContract("contract fingerprint mismatch")
        }
    }

    private static func fingerprint(_ contract: VivoSurrogateContract) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return SHA256.hash(data: try encoder.encode(contract)).map { String(format: "%02x", $0) }.joined()
    }
}

public struct VivoSurrogateBatch: Sendable, Equatable {
    public let batchSize: Int
    public let inputWidth: Int
    public let values: [Float]
    public init(batchSize: Int, inputWidth: Int, values: [Float]) throws {
        let expected = batchSize.multipliedReportingOverflow(by: inputWidth)
        guard batchSize > 0, inputWidth > 0, !expected.overflow,
              values.count == expected.partialValue, values.allSatisfy(\.isFinite) else {
            throw VivoSurrogateError.invalidBatch
        }
        self.batchSize = batchSize; self.inputWidth = inputWidth; self.values = values
    }
    public func fingerprint() -> String {
        var hasher = SHA256()
        for size in [batchSize, inputWidth] {
            var value = UInt64(size).littleEndian
            withUnsafeBytes(of: &value) { hasher.update(bufferPointer: $0) }
        }
        for element in values {
            var value = element.bitPattern.littleEndian
            withUnsafeBytes(of: &value) { hasher.update(bufferPointer: $0) }
        }
        return hasher.finalize().map { String(format: "%02x", $0) }.joined()
    }
}

public struct VivoSurrogatePrediction: Sendable, Equatable {
    public let batchSize: Int
    public let outputWidth: Int
    public let values: [Float]
    public let normalizedUncertainty: [Float]
    public let backend: String
    public let modelFingerprint: String
    public init(batchSize: Int, outputWidth: Int, values: [Float], normalizedUncertainty: [Float],
                backend: String, modelFingerprint: String) throws {
        let count = batchSize.multipliedReportingOverflow(by: outputWidth)
        guard batchSize > 0, outputWidth > 0, !count.overflow, values.count == count.partialValue,
              values.allSatisfy(\.isFinite), !backend.isEmpty, !modelFingerprint.isEmpty,
              normalizedUncertainty.count == batchSize || normalizedUncertainty.count == count.partialValue,
              normalizedUncertainty.allSatisfy({ $0.isFinite && $0 >= 0 }) else {
            throw VivoSurrogateError.invalidPrediction
        }
        self.batchSize = batchSize; self.outputWidth = outputWidth; self.values = values
        self.normalizedUncertainty = normalizedUncertainty; self.backend = backend; self.modelFingerprint = modelFingerprint
    }
}

public protocol VivoSurrogateBackend: Sendable {
    var backendID: String { get }
    var modelFingerprint: String { get }
    var inputWidth: Int { get }
    var outputWidth: Int { get }
    func predict(_ batch: VivoSurrogateBatch) async throws -> VivoSurrogatePrediction
}

/// The trusted authoritative participant must evaluate this exact request.
/// A request/receipt establishes execution provenance, not experimental validation.
public struct VivoSurrogateAuthorityRequest: Codable, Sendable, Equatable {
    public let identifier: UUID
    public let contractFingerprint: String
    public let modelFingerprint: String
    public let hostContextFingerprint: String?
    public let batchFingerprint: String
    public let batchSize: Int
    public let generation: UInt64
}

public enum VivoSurrogateDecisionReason: String, Codable, Sendable {
    case accepted, authorityInterval, modelFingerprintMismatch, inputShapeMismatch, inputOutOfDomain
    case excessiveExtrapolation, excessiveUncertainty, invalidOutput, outputOutOfBounds, backendFailure
    case unavailableUncertainty, evaluationInFlight
}

public struct VivoSurrogateDecision: Sendable, Equatable {
    public let accepted: Bool
    public let reason: VivoSurrogateDecisionReason
    public let prediction: VivoSurrogatePrediction?
    public let maximumObservedUncertainty: Float
    public let maximumObservedExtrapolation: Float
    public let requiresAuthoritativeEvaluation: Bool
    public let authorityRequest: VivoSurrogateAuthorityRequest?
}

/// Fail-closed, single-flight gate. New gates require an initial authoritative
/// evaluation. Pending refresh survives repeated calls and backend failures.
/// Inputs are in the declared physical feature units; normalization is a backend
/// responsibility. Uncertainty numbers must use the contract's declared metric.
public actor VivoSurrogateAuthorityGate {
    public let contract: VivoSurrogateContract
    private let backend: any VivoSurrogateBackend
    private var consecutiveAcceptedSteps: UInt32 = 0
    private var generation: UInt64 = 0
    private var needsRefresh = true
    private var inFlight = false
    private var pending: VivoSurrogateAuthorityRequest?
    public private(set) var lastAuthoritativeResultFingerprint: String?

    public init(contract: VivoSurrogateContract, backend: any VivoSurrogateBackend) throws {
        try contract.validateIntegrity()
        guard backend.modelFingerprint == contract.modelFingerprint, backend.inputWidth == contract.inputs.count,
              backend.outputWidth == contract.outputs.count else { throw VivoSurrogateError.backendContractMismatch }
        self.contract = contract; self.backend = backend
    }

    public func evaluate(_ batch: VivoSurrogateBatch) async -> VivoSurrogateDecision {
        if inFlight { return decision(.evaluationInFlight, requiresAuthority: false) }
        guard batch.inputWidth == contract.inputs.count else { return decision(.inputShapeMismatch) }
        if pending != nil { return decision(.authorityInterval) }
        if needsRefresh || consecutiveAcceptedSteps >= min(contract.mandatoryAuthorityInterval, contract.maximumConsecutiveAcceptedSteps) {
            requestAuthority(for: batch)
            return decision(.authorityInterval)
        }
        guard contract.uncertaintyKind != .none else {
            requestAuthority(for: batch)
            return decision(.unavailableUncertainty)
        }
        var extrapolation: Double = 0
        for row in 0..<batch.batchSize {
            for column in 0..<batch.inputWidth {
                let feature = contract.inputs[column]
                let value = Double(batch.values[row * batch.inputWidth + column])
                let lower = Double(feature.minimum), upper = Double(feature.maximum)
                let outside = max(0, lower - value, value - upper)
                let range = upper - lower
                if outside > 0 { extrapolation = max(extrapolation, range > 0 ? outside / range : Double.infinity) }
            }
        }
        let reportExtrapolation = Float(min(extrapolation, Double(Float.greatestFiniteMagnitude)))
        guard extrapolation <= Double(contract.maximumNormalizedExtrapolation) else {
            requestAuthority(for: batch)
            return decision(.excessiveExtrapolation, extrapolation: reportExtrapolation)
        }
        // Actors are reentrant at await. Reserve before calling any backend.
        inFlight = true
        defer { inFlight = false }
        let prediction: VivoSurrogatePrediction
        do { prediction = try await backend.predict(batch) }
        catch {
            requestAuthority(for: batch)
            return decision(.backendFailure, extrapolation: reportExtrapolation)
        }
        guard !Task.isCancelled, prediction.modelFingerprint == contract.modelFingerprint,
              prediction.backend == backend.backendID, prediction.batchSize == batch.batchSize,
              prediction.outputWidth == contract.outputs.count else {
            requestAuthority(for: batch)
            return decision(.invalidOutput, extrapolation: reportExtrapolation)
        }
        let uncertainty = prediction.normalizedUncertainty.max() ?? Float.greatestFiniteMagnitude
        guard uncertainty <= contract.maximumNormalizedUncertainty else {
            requestAuthority(for: batch)
            return decision(.excessiveUncertainty, uncertainty: uncertainty, extrapolation: reportExtrapolation)
        }
        guard outputIsValid(prediction.values, batchSize: batch.batchSize) else {
            requestAuthority(for: batch)
            return decision(.outputOutOfBounds, uncertainty: uncertainty, extrapolation: reportExtrapolation)
        }
        consecutiveAcceptedSteps += 1
        return .init(accepted: true, reason: .accepted, prediction: prediction,
                     maximumObservedUncertainty: uncertainty, maximumObservedExtrapolation: reportExtrapolation,
                     requiresAuthoritativeEvaluation: false, authorityRequest: nil)
    }

    /// Call only after a successful authoritative computation by the owning
    /// participant. Old/replayed/foreign receipts, invalid outputs and in-flight
    /// acknowledgements leave the pending request intact.
    public func recordAuthoritativeEvaluation(request: VivoSurrogateAuthorityRequest,
                                             resultFingerprint: String, values: [Float],
                                             numericalChecksPassed: Bool) throws {
        guard !inFlight, let pending, request == pending, numericalChecksPassed,
              resultFingerprint.count == 64,
              resultFingerprint.utf8.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }),
              outputIsValid(values, batchSize: request.batchSize), generation < UInt64.max else {
            throw VivoSurrogateError.invalidAuthorityReceipt
        }
        self.pending = nil; needsRefresh = false; consecutiveAcceptedSteps = 0; generation += 1
        lastAuthoritativeResultFingerprint = resultFingerprint
    }

    @available(*, unavailable, message: "Supply the pending request, authoritative result identity, values and numerical status.")
    public func recordAuthoritativeEvaluation() {}

    public func pendingAuthorityRequest() -> VivoSurrogateAuthorityRequest? { pending }

    private func requestAuthority(for batch: VivoSurrogateBatch) {
        needsRefresh = true
        if pending == nil {
            pending = .init(identifier: UUID(), contractFingerprint: contract.fingerprint,
                            modelFingerprint: contract.modelFingerprint, hostContextFingerprint: contract.hostContextFingerprint,
                            batchFingerprint: batch.fingerprint(), batchSize: batch.batchSize, generation: generation)
        }
    }

    private func outputIsValid(_ values: [Float], batchSize: Int) -> Bool {
        let count = batchSize.multipliedReportingOverflow(by: contract.outputs.count)
        guard !count.overflow, values.count == count.partialValue else { return false }
        for index in values.indices {
            let feature = contract.outputs[index % contract.outputs.count]
            guard values[index].isFinite, values[index] >= feature.minimum, values[index] <= feature.maximum else { return false }
        }
        return true
    }

    private func decision(_ reason: VivoSurrogateDecisionReason, uncertainty: Float = 0,
                          extrapolation: Float = 0, requiresAuthority: Bool = true) -> VivoSurrogateDecision {
        .init(accepted: false, reason: reason, prediction: nil, maximumObservedUncertainty: uncertainty,
              maximumObservedExtrapolation: extrapolation, requiresAuthoritativeEvaluation: requiresAuthority,
              authorityRequest: pending)
    }
}

public enum VivoSurrogateError: Error, LocalizedError, Sendable {
    case invalidContract(String), invalidBatch, invalidPrediction, backendContractMismatch
    case backendUnavailable(String), modelFailure(String), invalidAuthorityReceipt
    public var errorDescription: String? {
        switch self {
        case .invalidContract(let reason): return "invalid surrogate contract: \(reason)"
        case .invalidBatch: return "invalid surrogate input batch"
        case .invalidPrediction: return "invalid surrogate prediction"
        case .backendContractMismatch: return "surrogate backend does not match its contract"
        case .backendUnavailable(let reason): return "surrogate backend unavailable: \(reason)"
        case .modelFailure(let reason): return "surrogate model failed: \(reason)"
        case .invalidAuthorityReceipt: return "authoritative result is absent, invalid, stale or belongs to another request"
        }
    }
}

public struct VivoAffineSurrogateBackend: VivoSurrogateBackend {
    public let backendID = "numivivo.affine.reference.v1"
    public let modelFingerprint: String
    public let inputWidth: Int
    public let outputWidth: Int
    public let weights: [Float]
    public let bias: [Float]
    public let outputUncertainty: [Float]
    public init(modelFingerprint: String, inputWidth: Int, outputWidth: Int, weights: [Float],
                bias: [Float], outputUncertainty: [Float]) throws {
        let count = inputWidth.multipliedReportingOverflow(by: outputWidth)
        guard !modelFingerprint.isEmpty, inputWidth > 0, outputWidth > 0, !count.overflow,
              weights.count == count.partialValue, bias.count == outputWidth, outputUncertainty.count == outputWidth,
              weights.allSatisfy(\.isFinite), bias.allSatisfy(\.isFinite),
              outputUncertainty.allSatisfy({ $0.isFinite && $0 >= 0 }) else { throw VivoSurrogateError.backendContractMismatch }
        self.modelFingerprint = modelFingerprint; self.inputWidth = inputWidth; self.outputWidth = outputWidth
        self.weights = weights; self.bias = bias; self.outputUncertainty = outputUncertainty
    }
    public func predict(_ batch: VivoSurrogateBatch) async throws -> VivoSurrogatePrediction {
        let count = batch.batchSize.multipliedReportingOverflow(by: outputWidth)
        guard batch.inputWidth == inputWidth, !count.overflow else { throw VivoSurrogateError.invalidBatch }
        var output = [Float](repeating: 0, count: count.partialValue)
        var uncertainty = output
        for row in 0..<batch.batchSize {
            for destination in 0..<outputWidth {
                var accumulator = bias[destination]
                for source in 0..<inputWidth {
                    accumulator += weights[destination * inputWidth + source] * batch.values[row * inputWidth + source]
                }
                output[row * outputWidth + destination] = accumulator
                uncertainty[row * outputWidth + destination] = outputUncertainty[destination]
            }
        }
        return try .init(batchSize: batch.batchSize, outputWidth: outputWidth, values: output,
                         normalizedUncertainty: uncertainty, backend: backendID, modelFingerprint: modelFingerprint)
    }
}
