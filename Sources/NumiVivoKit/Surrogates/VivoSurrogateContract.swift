import CryptoKit
import Foundation

public struct VivoSurrogateFeature: Codable, Hashable, Sendable {
    public let id: String
    public let unit: String
    public let minimum: Float
    public let maximum: Float
    public let normalizationOffset: Float
    public let normalizationScale: Float

    public init(
        id: String,
        unit: String,
        minimum: Float,
        maximum: Float,
        normalizationOffset: Float = 0,
        normalizationScale: Float = 1
    ) {
        self.id = id
        self.unit = unit
        self.minimum = minimum
        self.maximum = maximum
        self.normalizationOffset = normalizationOffset
        self.normalizationScale = normalizationScale
    }

    public func normalize(_ value: Float) -> Float {
        (value - normalizationOffset) * normalizationScale
    }

    public func denormalize(_ value: Float) -> Float {
        value / normalizationScale + normalizationOffset
    }
}

public enum VivoSurrogateUncertaintyKind: String, Codable, Sendable {
    case none
    case predictedStandardDeviation
    case ensembleVariance
    case conformalRadius
}

/// Surrogate execution is permitted only inside this declared domain. The
/// contract never grants authority to bypass a numerical or biological monitor.
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

    public init(
        id: String,
        modelFingerprint: String,
        trainingDataFingerprint: String,
        mechanismPackFingerprint: String,
        hostContextFingerprint: String? = nil,
        inputs: [VivoSurrogateFeature],
        outputs: [VivoSurrogateFeature],
        uncertaintyKind: VivoSurrogateUncertaintyKind,
        maximumNormalizedUncertainty: Float,
        maximumNormalizedExtrapolation: Float = 0,
        maximumConsecutiveAcceptedSteps: UInt32 = 32,
        mandatoryAuthorityInterval: UInt32 = 32
    ) throws {
        let unsigned = VivoSurrogateContract(
            schemaVersion: 1,
            id: id,
            modelFingerprint: modelFingerprint,
            trainingDataFingerprint: trainingDataFingerprint,
            mechanismPackFingerprint: mechanismPackFingerprint,
            hostContextFingerprint: hostContextFingerprint,
            inputs: inputs,
            outputs: outputs,
            uncertaintyKind: uncertaintyKind,
            maximumNormalizedUncertainty: maximumNormalizedUncertainty,
            maximumNormalizedExtrapolation: maximumNormalizedExtrapolation,
            maximumConsecutiveAcceptedSteps: maximumConsecutiveAcceptedSteps,
            mandatoryAuthorityInterval: mandatoryAuthorityInterval,
            fingerprint: ""
        )
        try unsigned.validate()
        self = .init(
            schemaVersion: unsigned.schemaVersion,
            id: unsigned.id,
            modelFingerprint: unsigned.modelFingerprint,
            trainingDataFingerprint: unsigned.trainingDataFingerprint,
            mechanismPackFingerprint: unsigned.mechanismPackFingerprint,
            hostContextFingerprint: unsigned.hostContextFingerprint,
            inputs: unsigned.inputs,
            outputs: unsigned.outputs,
            uncertaintyKind: unsigned.uncertaintyKind,
            maximumNormalizedUncertainty: unsigned.maximumNormalizedUncertainty,
            maximumNormalizedExtrapolation: unsigned.maximumNormalizedExtrapolation,
            maximumConsecutiveAcceptedSteps: unsigned.maximumConsecutiveAcceptedSteps,
            mandatoryAuthorityInterval: unsigned.mandatoryAuthorityInterval,
            fingerprint: try Self.fingerprint(unsigned)
        )
    }

    private init(
        schemaVersion: UInt32,
        id: String,
        modelFingerprint: String,
        trainingDataFingerprint: String,
        mechanismPackFingerprint: String,
        hostContextFingerprint: String?,
        inputs: [VivoSurrogateFeature],
        outputs: [VivoSurrogateFeature],
        uncertaintyKind: VivoSurrogateUncertaintyKind,
        maximumNormalizedUncertainty: Float,
        maximumNormalizedExtrapolation: Float,
        maximumConsecutiveAcceptedSteps: UInt32,
        mandatoryAuthorityInterval: UInt32,
        fingerprint: String
    ) {
        self.schemaVersion = schemaVersion
        self.id = id
        self.modelFingerprint = modelFingerprint
        self.trainingDataFingerprint = trainingDataFingerprint
        self.mechanismPackFingerprint = mechanismPackFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.inputs = inputs
        self.outputs = outputs
        self.uncertaintyKind = uncertaintyKind
        self.maximumNormalizedUncertainty = maximumNormalizedUncertainty
        self.maximumNormalizedExtrapolation = maximumNormalizedExtrapolation
        self.maximumConsecutiveAcceptedSteps = maximumConsecutiveAcceptedSteps
        self.mandatoryAuthorityInterval = mandatoryAuthorityInterval
        self.fingerprint = fingerprint
    }

    public func validate() throws {
        guard schemaVersion == 1,
              !id.isEmpty,
              !modelFingerprint.isEmpty,
              !trainingDataFingerprint.isEmpty,
              !mechanismPackFingerprint.isEmpty,
              !inputs.isEmpty,
              !outputs.isEmpty,
              Set(inputs.map(\.id)).count == inputs.count,
              Set(outputs.map(\.id)).count == outputs.count,
              maximumNormalizedUncertainty.isFinite,
              maximumNormalizedUncertainty >= 0,
              maximumNormalizedExtrapolation.isFinite,
              maximumNormalizedExtrapolation >= 0,
              maximumConsecutiveAcceptedSteps > 0,
              mandatoryAuthorityInterval > 0,
              mandatoryAuthorityInterval <= maximumConsecutiveAcceptedSteps else {
            throw VivoSurrogateError.invalidContract("top-level contract is inconsistent")
        }
        for feature in inputs + outputs {
            guard !feature.id.isEmpty,
                  !feature.unit.isEmpty,
                  feature.minimum.isFinite,
                  feature.maximum.isFinite,
                  feature.minimum <= feature.maximum,
                  feature.normalizationOffset.isFinite,
                  feature.normalizationScale.isFinite,
                  feature.normalizationScale != 0 else {
                throw VivoSurrogateError.invalidContract("feature \(feature.id) has invalid bounds or normalization")
            }
        }
        if uncertaintyKind == .none, maximumNormalizedUncertainty != 0 {
            throw VivoSurrogateError.invalidContract("uncertainty threshold must be zero when uncertaintyKind is none")
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
        guard batchSize > 0,
              inputWidth > 0,
              !expected.overflow,
              values.count == expected.partialValue,
              values.allSatisfy(\.isFinite) else {
            throw VivoSurrogateError.invalidBatch
        }
        self.batchSize = batchSize
        self.inputWidth = inputWidth
        self.values = values
    }
}

public struct VivoSurrogatePrediction: Sendable, Equatable {
    public let batchSize: Int
    public let outputWidth: Int
    public let values: [Float]
    public let normalizedUncertainty: [Float]
    public let backend: String
    public let modelFingerprint: String

    public init(
        batchSize: Int,
        outputWidth: Int,
        values: [Float],
        normalizedUncertainty: [Float],
        backend: String,
        modelFingerprint: String
    ) throws {
        let outputElements = batchSize.multipliedReportingOverflow(by: outputWidth)
        guard batchSize > 0,
              outputWidth > 0,
              !outputElements.overflow,
              values.count == outputElements.partialValue,
              values.allSatisfy(\.isFinite),
              (normalizedUncertainty.count == batchSize || normalizedUncertainty.count == outputElements.partialValue),
              normalizedUncertainty.allSatisfy({ $0.isFinite && $0 >= 0 }) else {
            throw VivoSurrogateError.invalidPrediction
        }
        self.batchSize = batchSize
        self.outputWidth = outputWidth
        self.values = values
        self.normalizedUncertainty = normalizedUncertainty
        self.backend = backend
        self.modelFingerprint = modelFingerprint
    }
}

public protocol VivoSurrogateBackend: Sendable {
    var backendID: String { get }
    var modelFingerprint: String { get }
    var inputWidth: Int { get }
    var outputWidth: Int { get }

    func predict(_ batch: VivoSurrogateBatch) async throws -> VivoSurrogatePrediction
}

public enum VivoSurrogateDecisionReason: String, Codable, Sendable {
    case accepted
    case authorityInterval
    case modelFingerprintMismatch
    case inputShapeMismatch
    case inputOutOfDomain
    case excessiveExtrapolation
    case excessiveUncertainty
    case invalidOutput
    case outputOutOfBounds
    case backendFailure
}

public struct VivoSurrogateDecision: Sendable, Equatable {
    public let accepted: Bool
    public let reason: VivoSurrogateDecisionReason
    public let prediction: VivoSurrogatePrediction?
    public let maximumObservedUncertainty: Float
    public let maximumObservedExtrapolation: Float
    public let requiresAuthoritativeEvaluation: Bool
}

/// Central acceptance gate. It validates both input and output domains, enforces
/// mandatory authoritative refresh, and never converts a surrogate prediction
/// into a safety decision.
public actor VivoSurrogateAuthorityGate {
    public let contract: VivoSurrogateContract
    private let backend: any VivoSurrogateBackend
    private var consecutiveAcceptedSteps: UInt32 = 0

    public init(contract: VivoSurrogateContract, backend: any VivoSurrogateBackend) throws {
        try contract.validate()
        guard backend.modelFingerprint == contract.modelFingerprint,
              backend.inputWidth == contract.inputs.count,
              backend.outputWidth == contract.outputs.count else {
            throw VivoSurrogateError.backendContractMismatch
        }
        self.contract = contract
        self.backend = backend
    }

    public func evaluate(_ batch: VivoSurrogateBatch) async -> VivoSurrogateDecision {
        guard batch.inputWidth == contract.inputs.count else {
            return reject(.inputShapeMismatch)
        }
        if consecutiveAcceptedSteps >= contract.mandatoryAuthorityInterval {
            consecutiveAcceptedSteps = 0
            return reject(.authorityInterval)
        }

        var maximumExtrapolation: Float = 0
        for row in 0..<batch.batchSize {
            for column in 0..<batch.inputWidth {
                let feature = contract.inputs[column]
                let value = batch.values[row * batch.inputWidth + column]
                let range = max(feature.maximum - feature.minimum, Float.leastNonzeroMagnitude)
                let below = max(0, feature.minimum - value) / range
                let above = max(0, value - feature.maximum) / range
                maximumExtrapolation = max(maximumExtrapolation, below, above)
            }
        }
        guard maximumExtrapolation <= contract.maximumNormalizedExtrapolation else {
            return reject(.excessiveExtrapolation, extrapolation: maximumExtrapolation)
        }

        let prediction: VivoSurrogatePrediction
        do {
            prediction = try await backend.predict(batch)
        } catch {
            return reject(.backendFailure, extrapolation: maximumExtrapolation)
        }
        guard prediction.modelFingerprint == contract.modelFingerprint,
              prediction.batchSize == batch.batchSize,
              prediction.outputWidth == contract.outputs.count else {
            return reject(.invalidOutput, extrapolation: maximumExtrapolation)
        }

        let maximumUncertainty = prediction.normalizedUncertainty.max() ?? 0
        guard maximumUncertainty <= contract.maximumNormalizedUncertainty else {
            return reject(
                .excessiveUncertainty,
                uncertainty: maximumUncertainty,
                extrapolation: maximumExtrapolation
            )
        }
        for row in 0..<prediction.batchSize {
            for column in 0..<prediction.outputWidth {
                let value = prediction.values[row * prediction.outputWidth + column]
                let feature = contract.outputs[column]
                guard value.isFinite, value >= feature.minimum, value <= feature.maximum else {
                    return reject(
                        .outputOutOfBounds,
                        uncertainty: maximumUncertainty,
                        extrapolation: maximumExtrapolation
                    )
                }
            }
        }

        consecutiveAcceptedSteps &+= 1
        return .init(
            accepted: true,
            reason: .accepted,
            prediction: prediction,
            maximumObservedUncertainty: maximumUncertainty,
            maximumObservedExtrapolation: maximumExtrapolation,
            requiresAuthoritativeEvaluation: false
        )
    }

    public func recordAuthoritativeEvaluation() {
        consecutiveAcceptedSteps = 0
    }

    private func reject(
        _ reason: VivoSurrogateDecisionReason,
        uncertainty: Float = 0,
        extrapolation: Float = 0
    ) -> VivoSurrogateDecision {
        .init(
            accepted: false,
            reason: reason,
            prediction: nil,
            maximumObservedUncertainty: uncertainty,
            maximumObservedExtrapolation: extrapolation,
            requiresAuthoritativeEvaluation: true
        )
    }
}

public enum VivoSurrogateError: Error, LocalizedError, Sendable {
    case invalidContract(String)
    case invalidBatch
    case invalidPrediction
    case backendContractMismatch
    case backendUnavailable(String)
    case modelFailure(String)

    public var errorDescription: String? {
        switch self {
        case .invalidContract(let reason): return "invalid surrogate contract: \(reason)"
        case .invalidBatch: return "invalid surrogate input batch"
        case .invalidPrediction: return "invalid surrogate prediction"
        case .backendContractMismatch: return "surrogate backend does not match its contract"
        case .backendUnavailable(let reason): return "surrogate backend unavailable: \(reason)"
        case .modelFailure(let reason): return "surrogate model failed: \(reason)"
        }
    }
}

/// Deterministic reference backend for affine or locally linear surrogate blocks.
/// It is useful for compiler-generated reduced models and as a reference path for
/// validating Metal, Core ML, or external MLX implementations.
public struct VivoAffineSurrogateBackend: VivoSurrogateBackend {
    public let backendID = "numivivo.affine.reference.v1"
    public let modelFingerprint: String
    public let inputWidth: Int
    public let outputWidth: Int
    public let weights: [Float]
    public let bias: [Float]
    public let outputUncertainty: [Float]

    public init(
        modelFingerprint: String,
        inputWidth: Int,
        outputWidth: Int,
        weights: [Float],
        bias: [Float],
        outputUncertainty: [Float]
    ) throws {
        let matrixElements = inputWidth.multipliedReportingOverflow(by: outputWidth)
        guard inputWidth > 0,
              outputWidth > 0,
              !matrixElements.overflow,
              weights.count == matrixElements.partialValue,
              bias.count == outputWidth,
              outputUncertainty.count == outputWidth,
              weights.allSatisfy(\.isFinite),
              bias.allSatisfy(\.isFinite),
              outputUncertainty.allSatisfy({ $0.isFinite && $0 >= 0 }) else {
            throw VivoSurrogateError.backendContractMismatch
        }
        self.modelFingerprint = modelFingerprint
        self.inputWidth = inputWidth
        self.outputWidth = outputWidth
        self.weights = weights
        self.bias = bias
        self.outputUncertainty = outputUncertainty
    }

    public func predict(_ batch: VivoSurrogateBatch) async throws -> VivoSurrogatePrediction {
        guard batch.inputWidth == inputWidth else { throw VivoSurrogateError.invalidBatch }
        var output = [Float](repeating: 0, count: batch.batchSize * outputWidth)
        var uncertainty = [Float](repeating: 0, count: batch.batchSize * outputWidth)
        for row in 0..<batch.batchSize {
            for destination in 0..<outputWidth {
                var accumulator = bias[destination]
                let weightBase = destination * inputWidth
                let inputBase = row * inputWidth
                for source in 0..<inputWidth {
                    accumulator += weights[weightBase + source] * batch.values[inputBase + source]
                }
                output[row * outputWidth + destination] = accumulator
                uncertainty[row * outputWidth + destination] = outputUncertainty[destination]
            }
        }
        return try .init(
            batchSize: batch.batchSize,
            outputWidth: outputWidth,
            values: output,
            normalizedUncertainty: uncertainty,
            backend: backendID,
            modelFingerprint: modelFingerprint
        )
    }
}
