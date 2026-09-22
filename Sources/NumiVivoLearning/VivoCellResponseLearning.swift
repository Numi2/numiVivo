import Foundation
import MLX
import MLXNN
import MLXOptimizers
import NumiVivoKit

/// Errors raised by the MLX response learner before an artifact is published.
public enum VivoCellResponseLearningError: Error, LocalizedError, Sendable {
    case invalid(String)
    case incompatible(String)
    case limit(String)

    public var errorDescription: String? {
        switch self {
        case .invalid(let value): return "cell-response learning: \(value)"
        case .incompatible(let value): return "cell-response learning incompatible: \(value)"
        case .limit(let value): return "cell-response learning resource limit: \(value)"
        }
    }
}

/// The bounded neural architecture used for one source-bound response model.
/// Context cells, perturbation identity, and supplied target descriptors are
/// fused before separate mean and variance heads produce log-CPM outputs.
public struct VivoCellResponseArchitecture: Codable, Sendable, Equatable {
    public let hiddenWidth: Int
    public let contextCells: Int
    public let maximumParameters: Int

    public init(hiddenWidth: Int = 256, contextCells: Int = 32, maximumParameters: Int = 30_000_000) {
        self.hiddenWidth = hiddenWidth
        self.contextCells = contextCells
        self.maximumParameters = maximumParameters
    }

    private enum CodingKeys: String, CodingKey { case hiddenWidth, contextCells, maximumParameters }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["hiddenWidth", "contextCells", "maximumParameters"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        hiddenWidth = try values.decode(Int.self, forKey: .hiddenWidth)
        contextCells = try values.decode(Int.self, forKey: .contextCells)
        maximumParameters = try values.decode(Int.self, forKey: .maximumParameters)
    }

    private static func product(_ left: Int, _ right: Int) throws -> Int {
        let result = left.multipliedReportingOverflow(by: right)
        guard !result.overflow, result.partialValue >= 0 else {
            throw VivoCellResponseLearningError.limit("parameter arithmetic")
        }
        return result.partialValue
    }

    private static func sum(_ values: [Int]) throws -> Int {
        var total = 0
        for value in values {
            let result = total.addingReportingOverflow(value)
            guard !result.overflow else { throw VivoCellResponseLearningError.limit("parameter arithmetic") }
            total = result.partialValue
        }
        return total
    }

    func parameterCount(featureCount: Int, targetCount: Int, descriptorCount: Int) throws -> Int {
        let projection = try Self.product(featureCount, hiddenWidth)
        let embedding = try Self.product(targetCount, hiddenWidth)
        // A target-specific response prior is fitted exclusively from the
        // training partition. It gives known perturbations a measured starting
        // response instead of asking a randomly initialized full-gene decoder
        // to discover a 33k-dimensional effect from a few updates.
        let responsePrior = try Self.product(targetCount, featureCount)
        let descriptor = try Self.product(descriptorCount, hiddenWidth)
        let fusion = try Self.product(try Self.product(hiddenWidth, 3), hiddenWidth)
        let output = try Self.product(try Self.product(featureCount, hiddenWidth), 2)
        return try Self.sum([
            projection, hiddenWidth,
            embedding, responsePrior,
            descriptor, hiddenWidth,
            fusion, hiddenWidth,
            output, try Self.product(featureCount, 2)
        ])
    }

    func validate(featureCount: Int, targetCount: Int, descriptorCount: Int) throws {
        guard (16...1_024).contains(hiddenWidth), (1...512).contains(contextCells),
              (100_000...30_000_000).contains(maximumParameters),
              featureCount > 0, targetCount > 0, descriptorCount > 0 else {
            throw VivoCellResponseLearningError.invalid("architecture dimensions or parameter budget")
        }
        let count = try parameterCount(featureCount: featureCount, targetCount: targetCount, descriptorCount: descriptorCount)
        guard count <= maximumParameters else {
            throw VivoCellResponseLearningError.limit("architecture has \(count) parameters, above \(maximumParameters)")
        }
    }
}

/// A deterministic training run. The learner uses momentum-free SGD so a
/// checkpoint can resume exactly from the stored model, absolute step, and
/// stateless sampler seed.
public struct VivoCellResponseTrainingPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let architecture: VivoCellResponseArchitecture
    public let steps: Int
    public let batchSize: Int
    public let learningRate: Double
    public let weightDecay: Double
    public let seed: UInt64
    public let validationEvery: Int
    public let validationExamples: Int

    public init(id: String, architecture: VivoCellResponseArchitecture = .init(), steps: Int = 1_000,
                batchSize: Int = 4, learningRate: Double = 1e-3, weightDecay: Double = 1e-5,
                seed: UInt64 = 0, validationEvery: Int = 100, validationExamples: Int = 32) {
        schemaVersion = 1
        self.id = id
        self.architecture = architecture
        self.steps = steps
        self.batchSize = batchSize
        self.learningRate = learningRate
        self.weightDecay = weightDecay
        self.seed = seed
        self.validationEvery = validationEvery
        self.validationExamples = validationExamples
    }

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, id, architecture, steps, batchSize, learningRate, weightDecay, seed, validationEvery, validationExamples
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "schemaVersion", "id", "architecture", "steps", "batchSize", "learningRate", "weightDecay", "seed", "validationEvery", "validationExamples"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        architecture = try values.decode(VivoCellResponseArchitecture.self, forKey: .architecture)
        steps = try values.decode(Int.self, forKey: .steps)
        batchSize = try values.decode(Int.self, forKey: .batchSize)
        learningRate = try values.decode(Double.self, forKey: .learningRate)
        weightDecay = try values.decode(Double.self, forKey: .weightDecay)
        seed = try values.decode(UInt64.self, forKey: .seed)
        validationEvery = try values.decode(Int.self, forKey: .validationEvery)
        validationExamples = try values.decode(Int.self, forKey: .validationExamples)
    }

    func validateStatic() throws {
        guard schemaVersion == 1, !id.isEmpty, id.utf8.count <= 1_024,
              (1...1_000_000).contains(steps), (1...128).contains(batchSize),
              learningRate.isFinite, (1e-8...1).contains(learningRate),
              weightDecay.isFinite, (0...1).contains(weightDecay),
              (1...steps).contains(validationEvery), (1...4_096).contains(validationExamples) else {
            throw VivoCellResponseLearningError.invalid("training plan fields")
        }
        try architecture.validate(featureCount: 1, targetCount: 1, descriptorCount: 1)
    }

    func validate(for corpus: VivoCellResponseCorpusReader) throws {
        try validateStatic()
        try architecture.validate(featureCount: corpus.featureCount, targetCount: corpus.targetCount,
                                  descriptorCount: corpus.descriptorCount)
        let denseValues = Double(batchSize) * Double(architecture.contextCells) * Double(corpus.featureCount)
        guard denseValues.isFinite, denseValues <= 134_217_728 else {
            throw VivoCellResponseLearningError.limit("dense context batch exceeds 512 MiB")
        }
        let trainingExamples = try corpus.examples(in: .training)
        let trainingStrata = Set(trainingExamples.map(\.stratumIndex)).count
        let scheduled = steps.multipliedReportingOverflow(by: batchSize)
        guard !scheduled.overflow, scheduled.partialValue >= trainingStrata else {
            throw VivoCellResponseLearningError.invalid("training plan cannot cover every response stratum")
        }
        _ = try corpus.trainingTargetBindings()
    }
}

/// A resume request can only add deterministic steps. Architecture, optimizer
/// settings, and sampler seed remain bound to the parent checkpoint.
public struct VivoCellResponseResumePlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let additionalSteps: Int

    public init(id: String, additionalSteps: Int) {
        schemaVersion = 1
        self.id = id
        self.additionalSteps = additionalSteps
    }

    private enum CodingKeys: String, CodingKey { case schemaVersion, id, additionalSteps }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "additionalSteps"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        additionalSteps = try values.decode(Int.self, forKey: .additionalSteps)
    }

    func validate() throws {
        guard schemaVersion == 1, !id.isEmpty, id.utf8.count <= 1_024,
              (1...1_000_000).contains(additionalSteps) else {
            throw VivoCellResponseLearningError.invalid("resume plan fields")
        }
    }
}

/// A prediction asks for a response distribution from declared control cells.
/// A target can use a trained target embedding only when its exact ID was in
/// the fitted corpus; novel targets use their supplied descriptor vector alone.
public struct VivoCellResponsePredictionPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let target: VivoCellResponseTarget
    public let contextSampleIDs: [String]
    public let useTrainedTargetEmbedding: Bool

    public init(id: String, target: VivoCellResponseTarget, contextSampleIDs: [String],
                useTrainedTargetEmbedding: Bool = false) {
        schemaVersion = 1
        self.id = id
        self.target = target
        self.contextSampleIDs = contextSampleIDs
        self.useTrainedTargetEmbedding = useTrainedTargetEmbedding
    }

    private enum CodingKeys: String, CodingKey { case schemaVersion, id, target, contextSampleIDs, useTrainedTargetEmbedding }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "target", "contextSampleIDs", "useTrainedTargetEmbedding"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        target = try values.decode(VivoCellResponseTarget.self, forKey: .target)
        contextSampleIDs = try values.decode([String].self, forKey: .contextSampleIDs)
        useTrainedTargetEmbedding = try values.decode(Bool.self, forKey: .useTrainedTargetEmbedding)
    }

    public func validate(expectedDescriptorCount: Int? = nil) throws {
        try target.validate(expectedDescriptorCount: expectedDescriptorCount)
        guard schemaVersion == 1, vivoOmicsID(id), !contextSampleIDs.isEmpty,
              contextSampleIDs.count <= 4_096, Set(contextSampleIDs).count == contextSampleIDs.count,
              contextSampleIDs.allSatisfy(vivoOmicsID), target.descriptors.allSatisfy({ Float($0).isFinite }) else {
            throw VivoCellResponseLearningError.invalid("prediction query")
        }
    }
}

public struct VivoCellResponseTrainingMetrics: Codable, Sendable, Equatable {
    public let step: UInt64
    public let trainNegativeLogLikelihood: Double
    public let validationNegativeLogLikelihood: Double?
    public let validationRMSE: Double?
    /// RMSE of the exact matched-control mean on the same frozen validation
    /// rows as `validationRMSE`. This is a gate, not a fitted metric.
    public let validationMatchedControlRMSE: Double?
}

/// State carried in the native checkpoint. All fields that affect sampling or
/// numerical continuation are explicit rather than inferred from a filename.
public struct VivoCellResponseModelState: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let architecture: VivoCellResponseArchitecture
    public let featureAxis: VivoCellResponseFeatureAxis
    public let targetIDs: [String]
    public let trainedTargetBindings: [VivoCellResponseTrainedTargetBinding]
    public let descriptorCount: Int
    public let corpus: VivoFingerprint
    public let trainingPlan: VivoFingerprint
    public let step: UInt64
    public let seed: UInt64
    public let samplerVersion: UInt32
    /// Algorithm used to compute the frozen training-only response prior.
    public let targetResponsePriorVersion: UInt32
    /// Deterministic optimizer routing used for the learned response decoder.
    public let optimizerVersion: UInt32
    public let learningRate: Double
    public let weightDecay: Double
    public let latestMetrics: VivoCellResponseTrainingMetrics

    public init(schemaVersion: Int, architecture: VivoCellResponseArchitecture,
                featureAxis: VivoCellResponseFeatureAxis, targetIDs: [String],
                trainedTargetBindings: [VivoCellResponseTrainedTargetBinding], descriptorCount: Int,
                corpus: VivoFingerprint, trainingPlan: VivoFingerprint, step: UInt64, seed: UInt64,
                samplerVersion: UInt32, targetResponsePriorVersion: UInt32,
                optimizerVersion: UInt32,
                learningRate: Double, weightDecay: Double,
                latestMetrics: VivoCellResponseTrainingMetrics) {
        self.schemaVersion = schemaVersion
        self.architecture = architecture
        self.featureAxis = featureAxis
        self.targetIDs = targetIDs
        self.trainedTargetBindings = trainedTargetBindings
        self.descriptorCount = descriptorCount
        self.corpus = corpus
        self.trainingPlan = trainingPlan
        self.step = step
        self.seed = seed
        self.samplerVersion = samplerVersion
        self.targetResponsePriorVersion = targetResponsePriorVersion
        self.optimizerVersion = optimizerVersion
        self.learningRate = learningRate
        self.weightDecay = weightDecay
        self.latestMetrics = latestMetrics
    }

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, architecture, featureAxis, targetIDs, trainedTargetBindings, descriptorCount,
             corpus, trainingPlan, step, seed, samplerVersion, targetResponsePriorVersion,
             optimizerVersion,
             learningRate, weightDecay, latestMetrics
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "schemaVersion", "architecture", "featureAxis", "targetIDs", "trainedTargetBindings", "descriptorCount",
            "corpus", "trainingPlan", "step", "seed", "samplerVersion", "targetResponsePriorVersion",
            "optimizerVersion",
            "learningRate", "weightDecay", "latestMetrics"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        architecture = try values.decode(VivoCellResponseArchitecture.self, forKey: .architecture)
        featureAxis = try values.decode(VivoCellResponseFeatureAxis.self, forKey: .featureAxis)
        targetIDs = try values.decode([String].self, forKey: .targetIDs)
        trainedTargetBindings = try values.decode([VivoCellResponseTrainedTargetBinding].self, forKey: .trainedTargetBindings)
        descriptorCount = try values.decode(Int.self, forKey: .descriptorCount)
        corpus = try values.decode(VivoFingerprint.self, forKey: .corpus)
        trainingPlan = try values.decode(VivoFingerprint.self, forKey: .trainingPlan)
        step = try values.decode(UInt64.self, forKey: .step)
        seed = try values.decode(UInt64.self, forKey: .seed)
        samplerVersion = try values.decode(UInt32.self, forKey: .samplerVersion)
        targetResponsePriorVersion = try values.decode(UInt32.self, forKey: .targetResponsePriorVersion)
        optimizerVersion = try values.decode(UInt32.self, forKey: .optimizerVersion)
        learningRate = try values.decode(Double.self, forKey: .learningRate)
        weightDecay = try values.decode(Double.self, forKey: .weightDecay)
        latestMetrics = try values.decode(VivoCellResponseTrainingMetrics.self, forKey: .latestMetrics)
    }
}

public struct VivoCellResponseModelReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let checkpoint: VivoFingerprint
    public let trainingPlan: VivoFingerprint
    public let corpus: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// Distributional prediction on the corpus feature axis. `featureMask` marks
/// axis features measured in the query source; values where it is zero are
/// structural absence rather than a measured zero-expression prediction.
public struct VivoCellResponsePrediction: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let id: String
    public let featureAxis: VivoCellResponseFeatureAxis
    public let featureMask: [Float]
    /// The deterministic aggregate of the declared matched control cells.
    /// Keeping it separate makes the learned perturbation response auditable.
    public let contextBaselineLogCPM: [Float]
    /// The training-only target response prior plus learned residual over
    /// `contextBaselineLogCPM` (zero for descriptor-only queries).
    public let meanDeltaLogCPM: [Float]
    /// Exact source rows used to form the published control baseline.
    public let contextSourceRows: [Int]
    public let meanLogCPM: [Float]
    public let varianceLogCPM: [Float]
}

public struct VivoCellResponsePredictionReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let model: VivoFingerprint
    public let corpus: VivoFingerprint
    public let query: VivoFingerprint
    public let prediction: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// One immutable held-out target and the exact control rows used to score it.
public struct VivoCellResponseEvaluationExample: Codable, Sendable, Equatable {
    public let targetRow: Int
    public let contextRows: [Int]
}

/// A bounded, replayable software evaluation. It records the chosen rows so
/// reported metrics cannot drift with source ordering or a later sampler.
public struct VivoCellResponseEvaluation: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let partition: VivoCellResponsePartition
    public let maximumExamples: Int
    public let samplerVersion: UInt32
    public let seed: UInt64
    public let examples: Int
    public let observedFeatures: Int
    public let negativeLogLikelihood: Double
    public let rmse: Double
    /// RMSE of the exact declared control mean used as the prediction
    /// baseline, scored on exactly the rows and features in `selection`.
    public let matchedControlRMSE: Double
    public let selection: [VivoCellResponseEvaluationExample]
}

public struct VivoCellResponseEvaluationReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let model: VivoFingerprint
    public let corpus: VivoFingerprint
    public let evaluation: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// Apple-native MLX implementation. The context branch is order-invariant over
/// declared control cells; it predicts a treatment residual over their mean,
/// so the model does not have to relearn the unperturbed expression profile.
final class VivoCellResponseMLXModel: Module {
    @ModuleInfo var cellProjection: Linear
    @ModuleInfo var targetEmbedding: Embedding
    @ModuleInfo var targetResponsePrior: Embedding
    @ModuleInfo var descriptorProjection: Linear
    @ModuleInfo var fusion: Linear
    @ModuleInfo var meanHead: Linear
    @ModuleInfo var varianceHead: Linear

    let featureCount: Int
    let hiddenWidth: Int

    init(featureCount: Int, targetCount: Int, descriptorCount: Int, architecture: VivoCellResponseArchitecture,
         targetResponsePriorValues: MLXArray? = nil) {
        self.featureCount = featureCount
        self.hiddenWidth = architecture.hiddenWidth
        cellProjection = Linear(featureCount, architecture.hiddenWidth)
        targetEmbedding = Embedding(embeddingCount: targetCount, dimensions: architecture.hiddenWidth)
        let prior = targetResponsePriorValues ?? MLXArray(
            Array(repeating: Float(0), count: targetCount * featureCount), [targetCount, featureCount])
        // Preserve the empirical training-only prior in the normal checkpoint
        // parameter tree, while excluding it from `valueAndGrad` and the
        // optimizer. Learned layers can only add a residual correction.
        let frozenPrior = Embedding(weight: prior)
        frozenPrior.freeze()
        targetResponsePrior = frozenPrior
        descriptorProjection = Linear(descriptorCount, architecture.hiddenWidth)
        fusion = Linear(architecture.hiddenWidth * 3, architecture.hiddenWidth)
        // Start from the training-only empirical response prior. The residual
        // decoder learns only corrections, so it cannot begin by degrading a
        // matched-control-plus-response prediction with random full-axis noise.
        meanHead = Linear(
            weight: MLXArray(Array(repeating: Float(0), count: featureCount * architecture.hiddenWidth),
                             [featureCount, architecture.hiddenWidth]),
            bias: MLXArray(Array(repeating: Float(0), count: featureCount)))
        varianceHead = Linear(architecture.hiddenWidth, featureCount)
    }

    func outputs(context: MLXArray, targetIDs: MLXArray, knownTargetMask: MLXArray,
                 descriptors: MLXArray) -> (mean: MLXArray, variance: MLXArray,
                                            baseline: MLXArray, delta: MLXArray) {
        let batch = context.shape[0]
        let contexts = context.shape[1]
        let projected = relu(cellProjection(context.reshaped([batch * contexts, featureCount])))
            .reshaped([batch, contexts, hiddenWidth])
        let pooled = projected.mean(axis: 1)
        let baseline = context.mean(axis: 1)
        let target = targetEmbedding(targetIDs) * knownTargetMask
        let prior = targetResponsePrior(targetIDs) * knownTargetMask
        let descriptor = relu(descriptorProjection(descriptors))
        let state = relu(fusion(concatenated([pooled, target, descriptor], axis: 1)))
        let delta = prior + meanHead(state)
        return (baseline + delta, softplus(varianceHead(state)) + 1e-4, baseline, delta)
    }
}

/// Versioned MLX response-learning surface. Artifacts remain source-bound by
/// the corpus receipt and use the shared native checkpoint codec for weights.
public enum VivoCellResponseLearning {
    public static let format = "numivivo-cell-response-learning/v2"
    public static let modelFormat = "numivivo-cell-response-model/v5"
    public static let predictionFormat = "numivivo-cell-response-prediction/v5"
    public static let evaluationFormat = "numivivo-cell-response-evaluation/v2"
    static let samplerVersion: UInt32 = 3
    static let targetResponsePriorVersion: UInt32 = 1
    static let optimizerVersion: UInt32 = 1
    static let optimizerFormat = "sgd-momentum-0-mean-head-output-axis-v1"
}

private struct VivoCellResponseBatch {
    let context: MLXArray
    let targetIDs: MLXArray
    let knownTargetMask: MLXArray
    let descriptors: MLXArray
    let targets: MLXArray
    let featureMask: MLXArray
    let targetValues: [Float]
    let contextRows: [Int]
}

private struct VivoCellResponseLoadedArtifact {
    let model: VivoCellResponseMLXModel
    let state: VivoCellResponseModelState
    let receipt: VivoCellResponseModelReceipt
    let trainingPlan: VivoCellResponseTrainingPlan
    let trainingPlanBytes: Data
}

private enum VivoCellResponseArtifactIO {
    static let checkpointName = "checkpoint.nvckpt"
    static let planName = "training-plan.json"
    static let receiptName = "receipt.json"
    static let weightsSection = "weights.safetensors"
    static let stateSection = "state.json"

    static func requireNew(_ destination: URL) throws {
        let manager = FileManager.default
        guard destination.isFileURL, destination.path.utf8.count <= 8_192,
              !destination.path.contains("\0"), !manager.fileExists(atPath: destination.path) else {
            throw VivoCellResponseLearningError.invalid("output directory already exists or path is invalid")
        }
        var directory = ObjCBool(false)
        guard manager.fileExists(atPath: destination.deletingLastPathComponent().path, isDirectory: &directory),
              directory.boolValue else {
            throw VivoCellResponseLearningError.invalid("output parent directory")
        }
    }

    static func staging(_ parent: URL, prefix: String) throws -> URL {
        let directory = parent.appendingPathComponent("." + prefix + "-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: false,
                                                attributes: [.posixPermissions: 0o700])
        return directory
    }

    static func read(_ root: VivoRootedDocumentStore, _ name: String, maximumBytes: Int) throws -> Data {
        do {
            return try root.readDocument(name, maximumBytes: maximumBytes)
        } catch {
            throw VivoCellResponseLearningError.invalid("artifact document \(name)")
        }
    }

    static func write(_ data: Data, to root: URL, name: String, maximumBytes: Int) throws {
        guard data.count <= maximumBytes else { throw VivoCellResponseLearningError.limit("artifact document \(name)") }
        try data.write(to: root.appendingPathComponent(name), options: .withoutOverwriting)
    }

    static func corpusFingerprint(_ corpus: VivoCellResponseCorpusReader) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(corpus.receipt))
    }

    static func modelFingerprint(_ receipt: VivoCellResponseModelReceipt) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(receipt))
    }

    static func checkpoint(model: VivoCellResponseMLXModel, state: VivoCellResponseModelState,
                           parentCheckpoint: VivoFingerprint?) throws -> Data {
        let arrays = Dictionary(uniqueKeysWithValues: model.parameters().flattened())
        let metadata = [
            "format": VivoCellResponseLearning.modelFormat,
            "optimizer": VivoCellResponseLearning.optimizerFormat,
            "sampler": "splitmix64-v1",
            "targetResponsePrior": "balanced-training-mean-v1",
            "featureCount": String(state.featureAxis.featureIDs.count),
            "step": String(state.step)
        ]
        let weights = try MLX.saveToData(arrays: arrays, metadata: metadata)
        let request = VivoCheckpointBuildRequest(
            runtime: .cellResponseMLX,
            artifactFingerprint: state.corpus.hex,
            sourceFingerprint: state.corpus.hex,
            experimentFingerprint: state.trainingPlan.hex,
            parentCheckpointFingerprint: parentCheckpoint?.hex,
            stepIndex: state.step,
            logicalTime: Double(state.step),
            randomStreamVersion: state.samplerVersion,
            sections: [
                .init(id: weightsSection, encoding: .rawBytes, elementCount: UInt64(weights.count),
                      elementStride: 1, data: weights),
                try .canonicalJSON(id: stateSection, value: state)
            ],
            metadata: [
                "format": VivoCellResponseLearning.modelFormat,
                "optimizer": VivoCellResponseLearning.optimizerFormat,
                "sampler": "splitmix64-v1",
                "targetResponsePrior": "balanced-training-mean-v1"
            ])
        return try VivoCheckpointCodec.encode(request)
    }

    static func validateState(_ state: VivoCellResponseModelState,
                              trainingPlan: VivoCellResponseTrainingPlan) throws {
        guard state.schemaVersion == 4, state.samplerVersion == VivoCellResponseLearning.samplerVersion,
              state.targetResponsePriorVersion == VivoCellResponseLearning.targetResponsePriorVersion,
              state.optimizerVersion == VivoCellResponseLearning.optimizerVersion,
              state.featureAxis.featureIDs.count > 0, state.featureAxis.featureIDs.count <= 200_000,
              state.targetIDs.count > 0, Set(state.targetIDs).count == state.targetIDs.count,
              state.targetIDs.allSatisfy(vivoOmicsID), !state.trainedTargetBindings.isEmpty,
              state.trainedTargetBindings.count <= state.targetIDs.count,
              Set(state.trainedTargetBindings.map(\.id)).count == state.trainedTargetBindings.count,
              state.trainedTargetBindings.allSatisfy({ state.targetIDs.contains($0.id) }),
              state.descriptorCount > 0, state.learningRate.isFinite, state.learningRate > 0,
              state.weightDecay.isFinite, state.weightDecay >= 0,
              state.latestMetrics.step == state.step,
              state.latestMetrics.trainNegativeLogLikelihood.isFinite,
              state.latestMetrics.validationNegativeLogLikelihood?.isFinite ?? true,
              state.latestMetrics.validationRMSE?.isFinite ?? true,
              state.latestMetrics.validationRMSE.map({ $0 >= 0 }) ?? true,
              state.latestMetrics.validationMatchedControlRMSE?.isFinite ?? true,
              state.latestMetrics.validationMatchedControlRMSE.map({ $0 >= 0 }) ?? true,
              state.architecture == trainingPlan.architecture,
              state.seed == trainingPlan.seed,
              state.learningRate == trainingPlan.learningRate,
              state.weightDecay == trainingPlan.weightDecay else {
            throw VivoCellResponseLearningError.invalid("checkpoint state")
        }
        try state.featureAxis.validate()
        try state.architecture.validate(featureCount: state.featureAxis.featureIDs.count,
                                        targetCount: state.targetIDs.count,
                                        descriptorCount: state.descriptorCount)
    }

    static func load(_ directory: URL, implementation: VivoFingerprint) throws -> VivoCellResponseLoadedArtifact {
        let root: VivoRootedDocumentStore
        do {
            root = try .init(rootURL: directory)
        } catch {
            throw VivoCellResponseLearningError.invalid("model artifact directory")
        }
        let receiptBytes = try read(root, receiptName, maximumBytes: 1_048_576)
        let receipt = try VivoCanonicalJSON.decode(VivoCellResponseModelReceipt.self, from: receiptBytes)
        guard receipt.schemaVersion == 1, receipt.format == VivoCellResponseLearning.modelFormat,
              receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == receiptBytes else {
            throw VivoCellResponseLearningError.invalid("model receipt")
        }
        let trainingPlanBytes = try read(root, planName, maximumBytes: 4_194_304)
        let trainingPlan = try VivoCanonicalJSON.decode(VivoCellResponseTrainingPlan.self, from: trainingPlanBytes)
        guard try VivoCanonicalJSON.encode(trainingPlan) == trainingPlanBytes else {
            throw VivoCellResponseLearningError.invalid("noncanonical training plan")
        }
        try trainingPlan.validateStatic()
        let checkpointBytes = try read(root, checkpointName, maximumBytes: 1_073_741_824)
        guard try VivoCanonicalJSON.fingerprint(trainingPlanBytes) == receipt.trainingPlan,
              try VivoCanonicalJSON.fingerprint(checkpointBytes) == receipt.checkpoint else {
            throw VivoCellResponseLearningError.invalid("model artifact fingerprints")
        }
        let decoded = try VivoCheckpointCodec.decode(checkpointBytes)
        let stateBytes = try decoded.section(id: stateSection)
        let weightsBytes = try decoded.section(id: weightsSection)
        let state = try VivoCanonicalJSON.decode(VivoCellResponseModelState.self, from: stateBytes)
        let sections = decoded.manifest.sections
        let weightsInfo = sections.first(where: { $0.id == weightsSection })
        let stateInfo = sections.first(where: { $0.id == stateSection })
        guard decoded.manifest.runtime == .cellResponseMLX,
              decoded.manifest.artifactFingerprint == receipt.corpus.hex,
              decoded.manifest.sourceFingerprint == receipt.corpus.hex,
              decoded.manifest.experimentFingerprint == receipt.trainingPlan.hex,
              decoded.manifest.logicalTime == Double(state.step),
              decoded.manifest.randomStreamVersion == VivoCellResponseLearning.samplerVersion,
              sections.count == 2, Set(sections.map(\.id)) == Set([weightsSection, stateSection]),
              let weightsInfo, weightsInfo.encoding == .rawBytes,
              weightsInfo.length == UInt64(weightsBytes.count), weightsInfo.elementCount == UInt64(weightsBytes.count),
              weightsInfo.elementStride == 1,
              let stateInfo, stateInfo.encoding == .canonicalJSON,
              stateInfo.length == UInt64(stateBytes.count), stateInfo.elementCount == 1,
              stateInfo.elementStride == UInt32(clamping: stateBytes.count),
              decoded.manifest.metadata == [
                "format": VivoCellResponseLearning.modelFormat,
                "optimizer": VivoCellResponseLearning.optimizerFormat,
                "sampler": "splitmix64-v1",
                "targetResponsePrior": "balanced-training-mean-v1"
              ] else {
            throw VivoCellResponseLearningError.invalid("cell-response checkpoint manifest")
        }
        guard try VivoCanonicalJSON.encode(state) == stateBytes,
              state.corpus == receipt.corpus, state.trainingPlan == receipt.trainingPlan,
              state.step == decoded.manifest.stepIndex else {
            throw VivoCellResponseLearningError.invalid("cell-response checkpoint state binding")
        }
        try validateState(state, trainingPlan: trainingPlan)
        let (arrays, metadata) = try MLX.loadArraysAndMetadata(data: weightsBytes)
        guard metadata["format"] == VivoCellResponseLearning.modelFormat,
              metadata["optimizer"] == VivoCellResponseLearning.optimizerFormat,
              metadata["sampler"] == "splitmix64-v1",
              metadata["targetResponsePrior"] == "balanced-training-mean-v1",
              metadata["featureCount"] == String(state.featureAxis.featureIDs.count),
              metadata["step"] == String(state.step) else {
            throw VivoCellResponseLearningError.invalid("cell-response weight metadata")
        }
        let model = VivoCellResponseMLXModel(featureCount: state.featureAxis.featureIDs.count,
                                             targetCount: state.targetIDs.count,
                                             descriptorCount: state.descriptorCount,
                                             architecture: state.architecture)
        try model.update(parameters: ModuleParameters.unflattened(arrays), verify: .all)
        eval(model)
        for (_, parameter) in model.parameters().flattened() {
            guard parameter.asArray(Float.self).allSatisfy(\.isFinite) else {
                throw VivoCellResponseLearningError.invalid("nonfinite model parameter")
            }
        }
        return .init(model: model, state: state, receipt: receipt, trainingPlan: trainingPlan,
                     trainingPlanBytes: trainingPlanBytes)
    }

    static func publish(model: VivoCellResponseMLXModel, state: VivoCellResponseModelState,
                        trainingPlanBytes: Data, implementation: VivoFingerprint,
                        parentCheckpoint: VivoFingerprint?, to destination: URL) throws -> VivoCellResponseModelReceipt {
        try requireNew(destination)
        let trainingPlan = try VivoCanonicalJSON.decode(VivoCellResponseTrainingPlan.self, from: trainingPlanBytes)
        guard try VivoCanonicalJSON.encode(trainingPlan) == trainingPlanBytes else {
            throw VivoCellResponseLearningError.invalid("noncanonical training plan")
        }
        try trainingPlan.validateStatic()
        let temporary = try staging(destination.deletingLastPathComponent(), prefix: "numivivo-cell-response-model")
        defer { try? FileManager.default.removeItem(at: temporary) }
        let checkpoint = try checkpoint(model: model, state: state, parentCheckpoint: parentCheckpoint)
        let receipt = try VivoCellResponseModelReceipt(
            schemaVersion: 1, format: VivoCellResponseLearning.modelFormat,
            checkpoint: VivoCanonicalJSON.fingerprint(checkpoint),
            trainingPlan: VivoCanonicalJSON.fingerprint(trainingPlanBytes), corpus: state.corpus,
            implementation: implementation)
        try write(trainingPlanBytes, to: temporary, name: planName, maximumBytes: 4_194_304)
        try write(checkpoint, to: temporary, name: checkpointName, maximumBytes: 1_073_741_824)
        try write(VivoCanonicalJSON.encode(receipt), to: temporary, name: receiptName, maximumBytes: 1_048_576)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }
}

extension VivoCellResponseLearning {
    private static func withGPUExecution<R>(_ body: () throws -> R) rethrows -> R {
        return try Device.withDefaultDevice(.gpu) {
            try Stream.withNewDefaultStream(device: .gpu) {
                try body()
            }
        }
    }

    private static func mix(seed: UInt64, step: UInt64, lane: UInt64, salt: UInt64) -> UInt64 {
        var value = seed &+ 0x9e3779b97f4a7c15 &* (step &+ 1)
        value = value &+ 0xbf58476d1ce4e5b9 &* (lane &+ 1)
        value = value &+ salt
        value = (value ^ (value >> 30)) &* 0xbf58476d1ce4e5b9
        value = (value ^ (value >> 27)) &* 0x94d049bb133111eb
        return value ^ (value >> 31)
    }

    private static func sampledIndex(count: Int, seed: UInt64, step: UInt64, lane: UInt64, salt: UInt64) throws -> Int {
        guard count > 0 else { throw VivoCellResponseLearningError.invalid("empty sample population") }
        return Int(mix(seed: seed, step: step, lane: lane, salt: salt) % UInt64(count))
    }

    /// Select a reproducible subset without replacement. Sorting the selected
    /// indices makes traversal independent of hash-table iteration order.
    private static func sampledIndicesWithoutReplacement(count: Int, take: Int, seed: UInt64,
                                                         salt: UInt64) throws -> [Int] {
        guard count > 0, take > 0, take <= count else {
            throw VivoCellResponseLearningError.invalid("evaluation sample population")
        }
        var selected = Set<Int>()
        selected.reserveCapacity(take)
        for position in 0..<take {
            let upper = count - take + position
            let candidate = Int(mix(seed: seed, step: UInt64(position), lane: 0, salt: salt) % UInt64(upper + 1))
            if !selected.insert(candidate).inserted {
                selected.insert(upper)
            }
        }
        guard selected.count == take else { throw VivoCellResponseLearningError.invalid("evaluation subset") }
        return selected.sorted()
    }

    /// Group each perturbed source sample with its declared matched-control
    /// context. Sampling a group before a cell avoids letting high-cell-count
    /// guides dominate training or held-out metrics.
    private static func stratifiedExamples(_ examples: [VivoCellResponseTrainingExample]) throws -> [[VivoCellResponseTrainingExample]] {
        let groups = Dictionary(grouping: examples, by: \.stratumIndex)
        guard !groups.isEmpty else { throw VivoCellResponseLearningError.invalid("empty response strata") }
        return try groups.keys.sorted().map { index in
            guard let rows = groups[index], !rows.isEmpty else {
                throw VivoCellResponseLearningError.invalid("empty response stratum")
            }
            return rows.sorted { $0.targetRow < $1.targetRow }
        }
    }

    /// Choose a source-order-independent control context shared by train,
    /// evaluation, and inference. A fixed hash ranks the immutable row IDs;
    /// each row is consumed once before cycling only when the model width
    /// exceeds the declared control population.
    private static func deterministicContextRows(_ rows: [Int], count: Int) throws -> [Int] {
        let candidates = rows.sorted { left, right in
            let leftRank = mix(seed: 0x0d4e3c2b1a987654, step: UInt64(left), lane: 0, salt: 0x6841f29d)
            let rightRank = mix(seed: 0x0d4e3c2b1a987654, step: UInt64(right), lane: 0, salt: 0x6841f29d)
            return leftRank == rightRank ? left < right : leftRank < rightRank
        }
        guard count > 0, !candidates.isEmpty else {
            throw VivoCellResponseLearningError.invalid("empty deterministic control context")
        }
        return (0..<count).map { candidates[$0 % candidates.count] }
    }

    /// A replay-bound cyclic window through a deterministic control ordering.
    /// Unlike the inference/evaluation context, the starting point changes for
    /// each optimizer lane so training eventually sees the declared control
    /// population instead of repeatedly fitting one fixed subset.
    private static func trainingContextRows(_ rows: [Int], count: Int, seed: UInt64,
                                            step: UInt64, lane: UInt64, stratum: Int) throws -> [Int] {
        let candidates = rows.sorted { left, right in
            let leftRank = mix(seed: 0x0d4e3c2b1a987654, step: UInt64(left), lane: 0, salt: 0x6841f29d)
            let rightRank = mix(seed: 0x0d4e3c2b1a987654, step: UInt64(right), lane: 0, salt: 0x6841f29d)
            return leftRank == rightRank ? left < right : leftRank < rightRank
        }
        guard count > 0, !candidates.isEmpty, stratum >= 0 else {
            throw VivoCellResponseLearningError.invalid("empty training control context")
        }
        let start = try sampledIndex(count: candidates.count, seed: seed, step: step, lane: lane,
                                     salt: 0x721a5be9 ^ UInt64(stratum))
        return (0..<count).map { candidates[(start + $0) % candidates.count] }
    }

    /// The empirical perturbation response used to initialize known targets.
    /// For every training target×context stratum, it computes the difference
    /// between the mean treated log-CPM vector and the mean of its declared
    /// matched controls, then gives each context equal weight. Validation and
    /// test rows are never opened here.
    private static func trainingOnlyTargetResponsePrior(corpus: VivoCellResponseCorpusReader,
                                                        trainedTargetIndices: Set<Int>) throws -> MLXArray {
        let strata = try stratifiedExamples(try corpus.examples(in: .training))
        let features = corpus.featureCount
        var prior = [Double](repeating: 0, count: corpus.targetCount * features)
        var targetContexts = [Int](repeating: 0, count: corpus.targetCount)
        var controlMeans: [String: [Double]] = [:]

        func mean(_ rows: [Int]) throws -> [Double] {
            guard !rows.isEmpty else { throw VivoCellResponseLearningError.invalid("empty response-prior rows") }
            var result = [Double](repeating: 0, count: features)
            for (rowOffset, row) in rows.enumerated() {
                if rowOffset % 32 == 0 { try Task.checkCancellation() }
                let values = try corpus.normalizedRow(row)
                guard values.count == features else { throw VivoCellResponseLearningError.invalid("response-prior row width") }
                for feature in 0..<features where corpus.featureMask[feature] > 0 {
                    result[feature] += Double(values[feature])
                }
            }
            let divisor = Double(rows.count)
            return result.map { $0 / divisor }
        }

        for candidates in strata {
            try Task.checkCancellation()
            guard let first = candidates.first,
                  trainedTargetIndices.contains(first.targetIndex),
                  candidates.allSatisfy({ $0.targetIndex == first.targetIndex && $0.contextRows == first.contextRows }) else {
                throw VivoCellResponseLearningError.invalid("response-prior stratum")
            }
            let contextKey = first.contextRows.map(String.init).joined(separator: ",")
            let controlMean: [Double]
            if let existing = controlMeans[contextKey] {
                controlMean = existing
            } else {
                let computed = try mean(first.contextRows)
                controlMeans[contextKey] = computed
                controlMean = computed
            }
            let treatedMean = try mean(candidates.map(\.targetRow))
            let offset = first.targetIndex * features
            for feature in 0..<features where corpus.featureMask[feature] > 0 {
                prior[offset + feature] += treatedMean[feature] - controlMean[feature]
            }
            targetContexts[first.targetIndex] += 1
        }
        for target in trainedTargetIndices {
            guard targetContexts[target] > 0 else {
                throw VivoCellResponseLearningError.invalid("response prior has an unobserved target")
            }
            let divisor = Double(targetContexts[target])
            let offset = target * features
            for feature in 0..<features where corpus.featureMask[feature] > 0 {
                prior[offset + feature] /= divisor
            }
        }
        guard prior.allSatisfy(\.isFinite) else { throw VivoCellResponseLearningError.invalid("nonfinite response prior") }
        return MLXArray(prior.map(Float.init), [corpus.targetCount, features])
    }

    private static func trainedTargetIndices(corpus: VivoCellResponseCorpusReader,
                                             bindings: [VivoCellResponseTrainedTargetBinding]) throws -> Set<Int> {
        guard !bindings.isEmpty, bindings.count <= corpus.targetCount,
              Set(bindings.map(\.id)).count == bindings.count else {
            throw VivoCellResponseLearningError.incompatible("trained target bindings")
        }
        let vocabulary = Dictionary(uniqueKeysWithValues: corpus.plan.targets.enumerated().map { ($0.element.id, $0.offset) })
        var indices = Set<Int>()
        for binding in bindings {
            guard let index = vocabulary[binding.id],
                  try corpus.plan.targets[index].fingerprint() == binding.descriptorFingerprint else {
                throw VivoCellResponseLearningError.incompatible("trained target binding differs from corpus")
            }
            indices.insert(index)
        }
        return indices
    }

    private static func validateCheckpointCorpus(_ state: VivoCellResponseModelState,
                                                 corpus: VivoCellResponseCorpusReader) throws -> Set<Int> {
        guard state.featureAxis == corpus.plan.featureAxis,
              state.targetIDs == corpus.plan.targets.map(\.id),
              state.descriptorCount == corpus.descriptorCount else {
            throw VivoCellResponseLearningError.incompatible("checkpoint corpus schema differs")
        }
        let actualBindings = try corpus.trainingTargetBindings()
        guard state.trainedTargetBindings == actualBindings else {
            throw VivoCellResponseLearningError.incompatible("checkpoint trained targets differ from corpus")
        }
        return try trainedTargetIndices(corpus: corpus, bindings: state.trainedTargetBindings)
    }

    private static func loss(model: VivoCellResponseMLXModel, context: MLXArray, targetIDs: MLXArray,
                             knownTargetMask: MLXArray, descriptors: MLXArray, targets: MLXArray,
                             featureMask: MLXArray) -> MLXArray {
        let output = model.outputs(context: context, targetIDs: targetIDs, knownTargetMask: knownTargetMask,
                                   descriptors: descriptors)
        // Learn the perturbation response after subtracting the explicit
        // matched-control baseline. Algebraically this preserves the absolute
        // Gaussian likelihood while making the learned head a delta decoder.
        let response = targets - output.baseline
        let residual = output.delta - response
        let gaussianConstant = Float(log(2.0 * Double.pi))
        let terms = (residual * residual / output.variance + log(output.variance) + gaussianConstant) * 0.5
        return sum(terms * featureMask) / (sum(featureMask) * Float(context.shape[0]))
    }

    private static func trainingBatch(corpus: VivoCellResponseCorpusReader,
                                      strata: [[VivoCellResponseTrainingExample]],
                                      plan: VivoCellResponseTrainingPlan, trainedTargetIndices: Set<Int>,
                                      step: UInt64) throws -> VivoCellResponseBatch {
        let features = corpus.featureCount
        let contexts = plan.architecture.contextCells
        let batches = plan.batchSize
        guard corpus.featureMask.count == features, corpus.featureMask.contains(where: { $0 > 0 }) else {
            throw VivoCellResponseLearningError.invalid("corpus has no measured prediction features")
        }
        var contextValues: [Float] = []
        var targetValues: [Float] = []
        var targetIDs: [Int32] = []
        var known: [Float] = []
        var descriptors: [Float] = []
        contextValues.reserveCapacity(batches * contexts * features)
        targetValues.reserveCapacity(batches * features)
        targetIDs.reserveCapacity(batches); known.reserveCapacity(batches); descriptors.reserveCapacity(batches * corpus.descriptorCount)
        for lane in 0..<batches {
            // This enumeration covers every response stratum before cycling;
            // cell selection within a stratum remains deterministic but varied.
            let completed = step.subtractingReportingOverflow(1)
            let scaled = completed.partialValue.multipliedReportingOverflow(by: UInt64(batches))
            let scheduled = scaled.partialValue.addingReportingOverflow(UInt64(lane))
            guard !completed.overflow, !scaled.overflow, !scheduled.overflow else {
                throw VivoCellResponseLearningError.limit("training schedule")
            }
            let stratumIndex = Int(scheduled.partialValue % UInt64(strata.count))
            let candidates = strata[stratumIndex]
            let exampleIndex = try sampledIndex(count: candidates.count, seed: plan.seed, step: step,
                                                lane: UInt64(lane), salt: 0x4c7912ae)
            let example = candidates[exampleIndex]
            guard trainedTargetIndices.contains(example.targetIndex),
                  let targetID = Int32(exactly: example.targetIndex), example.contextRows.count > 0 else {
                throw VivoCellResponseLearningError.invalid("training target or context rows")
            }
            for row in try trainingContextRows(example.contextRows, count: contexts, seed: plan.seed,
                                               step: step, lane: UInt64(lane), stratum: stratumIndex) {
                contextValues.append(contentsOf: try corpus.normalizedRow(row))
            }
            targetValues.append(contentsOf: try corpus.normalizedRow(example.targetRow))
            targetIDs.append(targetID); known.append(1)
            descriptors.append(contentsOf: corpus.plan.targets[example.targetIndex].descriptors.map(Float.init))
        }
        return .init(
            context: MLXArray(contextValues, [batches, contexts, features]),
            targetIDs: MLXArray(targetIDs, [batches]),
            knownTargetMask: MLXArray(known, [batches, 1]),
            descriptors: MLXArray(descriptors, [batches, corpus.descriptorCount]),
            targets: MLXArray(targetValues, [batches, features]),
            featureMask: MLXArray(corpus.featureMask, [1, features]),
            targetValues: targetValues, contextRows: [])
    }

    private static func evaluationBatch(corpus: VivoCellResponseCorpusReader,
                                        example: VivoCellResponseTrainingExample,
                                        architecture: VivoCellResponseArchitecture,
                                        trainedTargetIndices: Set<Int>) throws -> VivoCellResponseBatch {
        guard example.contextRows.count > 0 else {
            throw VivoCellResponseLearningError.invalid("evaluation target or context rows")
        }
        let known = trainedTargetIndices.contains(example.targetIndex)
        let targetID: Int32
        if known {
            guard let value = Int32(exactly: example.targetIndex) else {
                throw VivoCellResponseLearningError.invalid("evaluation target index")
            }
            targetID = value
        } else {
            targetID = 0
        }
        let selectedRows = try deterministicContextRows(example.contextRows, count: architecture.contextCells)
        var contextValues: [Float] = []
        for row in selectedRows {
            contextValues.append(contentsOf: try corpus.normalizedRow(row))
        }
        let targetValues = try corpus.normalizedRow(example.targetRow)
        return .init(
            context: MLXArray(contextValues, [1, architecture.contextCells, corpus.featureCount]),
            targetIDs: MLXArray([targetID], [1]), knownTargetMask: MLXArray([known ? Float(1) : Float(0)], [1, 1]),
            descriptors: MLXArray(corpus.plan.targets[example.targetIndex].descriptors.map(Float.init), [1, corpus.descriptorCount]),
            targets: MLXArray(targetValues, [1, corpus.featureCount]),
            featureMask: MLXArray(corpus.featureMask, [1, corpus.featureCount]), targetValues: targetValues,
            contextRows: selectedRows)
    }

    private static func evaluate(model: VivoCellResponseMLXModel, corpus: VivoCellResponseCorpusReader,
                                 examples: [VivoCellResponseTrainingExample], architecture: VivoCellResponseArchitecture,
                                 trainedTargetIndices: Set<Int>, maximumExamples: Int, seed: UInt64,
                                 salt: UInt64) throws -> (Double, Double, Double, Int, [VivoCellResponseEvaluationExample]) {
        guard (1...4_096).contains(maximumExamples) else { throw VivoCellResponseLearningError.invalid("evaluation example count") }
        let strata = try stratifiedExamples(examples)
        let take = min(maximumExamples, strata.count)
        guard take > 0 else { throw VivoCellResponseLearningError.invalid("evaluation partition") }
        var nll = 0.0, squaredError = 0.0, baselineSquaredError = 0.0, observed = 0
        var selection: [VivoCellResponseEvaluationExample] = []
        for index in try sampledIndicesWithoutReplacement(count: strata.count, take: take, seed: seed, salt: salt) {
            let candidates = strata[index]
            let row = try sampledIndex(count: candidates.count, seed: seed, step: UInt64(index), lane: 0, salt: salt ^ 0x38b64ea1)
            let example = candidates[row]
            let batch = try evaluationBatch(corpus: corpus, example: example, architecture: architecture,
                                            trainedTargetIndices: trainedTargetIndices)
            selection.append(.init(targetRow: example.targetRow, contextRows: batch.contextRows))
            let output = model.outputs(context: batch.context, targetIDs: batch.targetIDs,
                                       knownTargetMask: batch.knownTargetMask, descriptors: batch.descriptors)
            let means = output.mean.asArray(Float.self), variances = output.variance.asArray(Float.self),
                baselines = output.baseline.asArray(Float.self)
            guard means.count == corpus.featureCount, variances.count == corpus.featureCount,
                  baselines.count == corpus.featureCount else {
                throw VivoCellResponseLearningError.invalid("evaluation output shape")
            }
            for feature in 0..<corpus.featureCount where corpus.featureMask[feature] > 0 {
                let mean = Double(means[feature]), variance = Double(variances[feature]),
                    baseline = Double(baselines[feature]), target = Double(batch.targetValues[feature])
                guard mean.isFinite, variance.isFinite, variance > 0, baseline.isFinite, target.isFinite else {
                    throw VivoCellResponseLearningError.invalid("nonfinite evaluation output")
                }
                let residual = mean - target
                nll += 0.5 * (residual * residual / variance + log(variance) + log(2.0 * Double.pi))
                squaredError += residual * residual
                let baselineResidual = baseline - target
                baselineSquaredError += baselineResidual * baselineResidual
                observed += 1
            }
        }
        guard observed > 0 else { throw VivoCellResponseLearningError.invalid("evaluation has no measured features") }
        return (nll / Double(observed), sqrt(squaredError / Double(observed)),
                sqrt(baselineSquaredError / Double(observed)), observed, selection)
    }

    private static func runSteps(model: VivoCellResponseMLXModel, corpus: VivoCellResponseCorpusReader,
                                 plan: VivoCellResponseTrainingPlan, startingStep: UInt64,
                                 additionalSteps: Int,
                                 trainedTargetIndices: Set<Int>) throws -> VivoCellResponseTrainingMetrics {
        let examples = try corpus.examples(in: .training)
        let strata = try stratifiedExamples(examples)
        let validation = (try? corpus.examples(in: .validation)) ?? []
        let observedFeatureCount = corpus.featureMask.reduce(into: 0) { total, value in
            if value > 0 { total += 1 }
        }
        guard observedFeatureCount > 0 else {
            throw VivoCellResponseLearningError.invalid("training corpus has no measured features")
        }
        let outputAxisScale = Float(observedFeatureCount)
        let baseLearningRate = Float(plan.learningRate)
        let meanHeadLearningRate = baseLearningRate * outputAxisScale
        let meanHeadWeightDecay = Float(plan.weightDecay) / outputAxisScale
        guard outputAxisScale.isFinite, baseLearningRate.isFinite,
              meanHeadLearningRate.isFinite, meanHeadLearningRate > 0,
              meanHeadWeightDecay.isFinite, meanHeadWeightDecay >= 0 else {
            throw VivoCellResponseLearningError.limit("mean-head optimizer scale")
        }
        // The likelihood remains averaged over measured genes. Scale only the
        // independent full-axis mean decoder back to a per-feature update;
        // scaling shared or variance parameters would multiply their update by
        // the whole transcriptome width.
        let optimizer = MultiOptimizer(
            optimizers: [
                SGD(learningRate: meanHeadLearningRate, momentum: 0, weightDecay: meanHeadWeightDecay),
                SGD(learningRate: baseLearningRate, momentum: 0, weightDecay: Float(plan.weightDecay))
            ],
            filters: [{ key, _ in key == "meanHead.weight" || key == "meanHead.bias" }])
        let lossAndGrad = valueAndGrad(model: model) { model, arrays in
            [loss(model: model, context: arrays[0], targetIDs: arrays[1], knownTargetMask: arrays[2],
                  descriptors: arrays[3], targets: arrays[4], featureMask: arrays[5])]
        }
        var latest: VivoCellResponseTrainingMetrics?
        for offset in 0..<additionalSteps {
            try Task.checkCancellation()
            let increment = UInt64(offset + 1)
            let stepResult = startingStep.addingReportingOverflow(increment)
            guard !stepResult.overflow else { throw VivoCellResponseLearningError.limit("checkpoint step") }
            let step = stepResult.partialValue
            let batch = try trainingBatch(corpus: corpus, strata: strata, plan: plan,
                                          trainedTargetIndices: trainedTargetIndices, step: step)
            let (values, gradients) = lossAndGrad(model, [batch.context, batch.targetIDs, batch.knownTargetMask,
                                                           batch.descriptors, batch.targets, batch.featureMask])
            optimizer.update(model: model, gradients: gradients)
            eval(model, optimizer)
            guard let scalar = values[0].asArray(Float.self).first, scalar.isFinite else {
                throw VivoCellResponseLearningError.invalid("nonfinite training loss")
            }
            var validationNLL: Double?, validationRMSE: Double?, validationMatchedControlRMSE: Double?
            if step % UInt64(plan.validationEvery) == 0 || offset + 1 == additionalSteps, !validation.isEmpty {
                let result = try evaluate(model: model, corpus: corpus, examples: validation, architecture: plan.architecture,
                                          trainedTargetIndices: trainedTargetIndices, maximumExamples: plan.validationExamples,
                                          seed: plan.seed, salt: 0x6f3a5c29)
                validationNLL = result.0; validationRMSE = result.1; validationMatchedControlRMSE = result.2
            }
            latest = .init(step: step, trainNegativeLogLikelihood: Double(scalar),
                           validationNegativeLogLikelihood: validationNLL, validationRMSE: validationRMSE,
                           validationMatchedControlRMSE: validationMatchedControlRMSE)
        }
        guard let latest else { throw VivoCellResponseLearningError.invalid("training had no steps") }
        return latest
    }

    private static func inferenceBatch(corpus: VivoCellResponseCorpusReader, state: VivoCellResponseModelState,
                                       plan: VivoCellResponsePredictionPlan) throws -> VivoCellResponseBatch {
        try plan.validate(expectedDescriptorCount: state.descriptorCount)
        guard corpus.plan.featureAxis == state.featureAxis,
              corpus.featureMask.count == state.featureAxis.featureIDs.count,
              corpus.featureMask.contains(where: { $0 > 0 }) else {
            throw VivoCellResponseLearningError.invalid("prediction query")
        }
        var rows: [Int] = []
        for sampleID in plan.contextSampleIDs.sorted() {
            let assignment = try corpus.assignment(forSampleID: sampleID)
            guard assignment.role == .control else {
                throw VivoCellResponseLearningError.invalid("prediction context sample is not a declared control")
            }
            rows.append(contentsOf: try corpus.rows(forSampleID: sampleID))
        }
        guard !rows.isEmpty else { throw VivoCellResponseLearningError.invalid("prediction context has no cells") }
        let targetIndex: Int32
        let knownTarget: Float
        if plan.useTrainedTargetEmbedding {
            guard let binding = state.trainedTargetBindings.first(where: { $0.id == plan.target.id }),
                  try plan.target.fingerprint() == binding.descriptorFingerprint,
                  let index = state.targetIDs.firstIndex(of: plan.target.id), let value = Int32(exactly: index) else {
                throw VivoCellResponseLearningError.incompatible("prediction target has no trained embedding")
            }
            targetIndex = value; knownTarget = 1
        } else {
            targetIndex = 0; knownTarget = 0
        }
        let selectedRows = try deterministicContextRows(rows, count: state.architecture.contextCells)
        var contexts: [Float] = []
        contexts.reserveCapacity(state.architecture.contextCells * corpus.featureCount)
        for row in selectedRows {
            contexts.append(contentsOf: try corpus.normalizedRow(row))
        }
        let descriptors = plan.target.descriptors.map(Float.init)
        return .init(
            context: MLXArray(contexts, [1, state.architecture.contextCells, corpus.featureCount]),
            targetIDs: MLXArray([targetIndex], [1]), knownTargetMask: MLXArray([knownTarget], [1, 1]),
            descriptors: MLXArray(descriptors, [1, state.descriptorCount]),
            targets: MLXArray([Float](repeating: 0, count: corpus.featureCount), [1, corpus.featureCount]),
            featureMask: MLXArray(corpus.featureMask, [1, corpus.featureCount]), targetValues: [],
            contextRows: selectedRows)
    }

    private static func publishPrediction(_ prediction: VivoCellResponsePrediction,
                                          query: VivoCellResponsePredictionPlan,
                                          modelReceipt: VivoCellResponseModelReceipt,
                                          corpus: VivoCellResponseCorpusReader,
                                          implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponsePredictionReceipt {
        try VivoCellResponseArtifactIO.requireNew(destination)
        let queryBytes = try VivoCanonicalJSON.encode(query)
        let predictionBytes = try VivoCanonicalJSON.encode(prediction)
        let receipt = try VivoCellResponsePredictionReceipt(
            schemaVersion: 1, format: predictionFormat,
            model: VivoCellResponseArtifactIO.modelFingerprint(modelReceipt),
            corpus: VivoCellResponseArtifactIO.corpusFingerprint(corpus),
            query: VivoCanonicalJSON.fingerprint(queryBytes), prediction: VivoCanonicalJSON.fingerprint(predictionBytes),
            implementation: implementation)
        let temporary = try VivoCellResponseArtifactIO.staging(destination.deletingLastPathComponent(), prefix: "numivivo-cell-response-prediction")
        defer { try? FileManager.default.removeItem(at: temporary) }
        try VivoCellResponseArtifactIO.write(queryBytes, to: temporary, name: "query.json", maximumBytes: 4_194_304)
        try VivoCellResponseArtifactIO.write(predictionBytes, to: temporary, name: "prediction.json", maximumBytes: 268_435_456)
        try VivoCellResponseArtifactIO.write(VivoCanonicalJSON.encode(receipt), to: temporary, name: "receipt.json", maximumBytes: 1_048_576)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }

    private static func publishEvaluation(_ evaluation: VivoCellResponseEvaluation,
                                          modelReceipt: VivoCellResponseModelReceipt,
                                          corpus: VivoCellResponseCorpusReader,
                                          implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseEvaluationReceipt {
        try validateEvaluation(evaluation)
        try VivoCellResponseArtifactIO.requireNew(destination)
        let evaluationBytes = try VivoCanonicalJSON.encode(evaluation)
        let receipt = VivoCellResponseEvaluationReceipt(
            schemaVersion: 1, format: evaluationFormat,
            model: try VivoCellResponseArtifactIO.modelFingerprint(modelReceipt),
            corpus: try VivoCellResponseArtifactIO.corpusFingerprint(corpus),
            evaluation: try VivoCanonicalJSON.fingerprint(evaluationBytes), implementation: implementation)
        let temporary = try VivoCellResponseArtifactIO.staging(destination.deletingLastPathComponent(), prefix: "numivivo-cell-response-evaluation")
        defer { try? FileManager.default.removeItem(at: temporary) }
        try VivoCellResponseArtifactIO.write(evaluationBytes, to: temporary, name: "evaluation.json", maximumBytes: 134_217_728)
        try VivoCellResponseArtifactIO.write(try VivoCanonicalJSON.encode(receipt), to: temporary, name: "receipt.json", maximumBytes: 1_048_576)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }

    private static func validateEvaluation(_ evaluation: VivoCellResponseEvaluation) throws {
        guard evaluation.schemaVersion == 2, evaluation.format == evaluationFormat,
              (1...4_096).contains(evaluation.maximumExamples),
              evaluation.samplerVersion == samplerVersion,
              evaluation.examples == evaluation.selection.count, !evaluation.selection.isEmpty,
              evaluation.observedFeatures > 0, evaluation.negativeLogLikelihood.isFinite,
              evaluation.rmse.isFinite, evaluation.rmse >= 0,
              evaluation.matchedControlRMSE.isFinite, evaluation.matchedControlRMSE >= 0,
              Set(evaluation.selection.map(\.targetRow)).count == evaluation.selection.count,
              evaluation.selection.allSatisfy({ $0.targetRow >= 0 && (1...512).contains($0.contextRows.count) && $0.contextRows.allSatisfy({ $0 >= 0 }) }) else {
            throw VivoCellResponseLearningError.invalid("evaluation artifact")
        }
    }

    /// Require a verified response evaluation to beat its exact matched-control
    /// baseline. Raw evaluations remain publishable for diagnosis; callers use
    /// this explicit gate before treating one as a qualified engineering result.
    public static func requireMatchedControlImprovement(_ evaluation: VivoCellResponseEvaluation) throws {
        try validateEvaluation(evaluation)
        guard evaluation.rmse < evaluation.matchedControlRMSE else {
            throw VivoCellResponseLearningError.invalid(
                "evaluation does not improve the exact matched-control baseline")
        }
    }

    private static func validatePrediction(_ prediction: VivoCellResponsePrediction) throws {
        try prediction.featureAxis.validate()
        guard prediction.schemaVersion == 1, prediction.format == predictionFormat, vivoOmicsID(prediction.id),
              prediction.featureAxis.featureIDs.count == prediction.featureMask.count,
              prediction.featureMask.count == prediction.meanLogCPM.count,
              prediction.featureMask.count == prediction.contextBaselineLogCPM.count,
              prediction.featureMask.count == prediction.meanDeltaLogCPM.count,
              prediction.meanLogCPM.count == prediction.varianceLogCPM.count,
              (1...512).contains(prediction.contextSourceRows.count),
              prediction.contextSourceRows.allSatisfy({ $0 >= 0 }),
              prediction.featureMask.allSatisfy({ $0 == 0 || $0 == 1 }),
              prediction.featureMask.contains(where: { $0 > 0 }),
              prediction.contextBaselineLogCPM.allSatisfy(\.isFinite),
              prediction.meanDeltaLogCPM.allSatisfy(\.isFinite),
              prediction.meanLogCPM.allSatisfy(\.isFinite),
              prediction.varianceLogCPM.allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoCellResponseLearningError.invalid("prediction output")
        }
        guard zip(zip(prediction.meanLogCPM, prediction.contextBaselineLogCPM), prediction.meanDeltaLogCPM).allSatisfy({ values in
            let expected = values.0.1 + values.1
            let scale = max(Float(1), max(abs(values.0.0), abs(expected)))
            return abs(values.0.0 - expected) <= 1e-5 * scale
        }) else {
            throw VivoCellResponseLearningError.invalid("prediction baseline and residual disagree")
        }
    }

    private static func approximatelyEqual(_ expected: [Float], _ actual: [Float]) -> Bool {
        guard expected.count == actual.count else { return false }
        return zip(expected, actual).allSatisfy { values in
            let scale = max(Float(1), max(abs(values.0), abs(values.1)))
            return abs(values.0 - values.1) <= 1e-5 * scale
        }
    }

    /// Fit a model only from a verified corpus reader. The caller must open the
    /// corpus against its exact source count store before this function runs.
    public static func train(corpus: VivoCellResponseCorpusReader, plan: VivoCellResponseTrainingPlan,
                             implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseModelReceipt {
        return try withGPUExecution {
            try plan.validate(for: corpus)
            try VivoCellResponseArtifactIO.requireNew(destination)
            let planBytes = try VivoCanonicalJSON.encode(plan)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            let trainedTargetBindings = try corpus.trainingTargetBindings()
            let trainedTargetIndices = try trainedTargetIndices(corpus: corpus, bindings: trainedTargetBindings)
            MLXRandom.seed(plan.seed)
            let responsePrior = try trainingOnlyTargetResponsePrior(corpus: corpus,
                                                                     trainedTargetIndices: trainedTargetIndices)
            let model = VivoCellResponseMLXModel(featureCount: corpus.featureCount, targetCount: corpus.targetCount,
                                                 descriptorCount: corpus.descriptorCount, architecture: plan.architecture,
                                                 targetResponsePriorValues: responsePrior)
            let metrics = try runSteps(model: model, corpus: corpus, plan: plan, startingStep: 0,
                                       additionalSteps: plan.steps, trainedTargetIndices: trainedTargetIndices)
            let state = VivoCellResponseModelState(
                schemaVersion: 4, architecture: plan.architecture, featureAxis: corpus.plan.featureAxis,
                targetIDs: corpus.plan.targets.map(\.id), trainedTargetBindings: trainedTargetBindings,
                descriptorCount: corpus.descriptorCount,
                corpus: corpusFingerprint, trainingPlan: try VivoCanonicalJSON.fingerprint(planBytes),
                step: metrics.step, seed: plan.seed, samplerVersion: samplerVersion,
                targetResponsePriorVersion: targetResponsePriorVersion,
                optimizerVersion: optimizerVersion,
                learningRate: plan.learningRate, weightDecay: plan.weightDecay, latestMetrics: metrics)
            return try VivoCellResponseArtifactIO.publish(model: model, state: state, trainingPlanBytes: planBytes,
                                                           implementation: implementation, parentCheckpoint: nil, to: destination)
        }
    }

    /// Continue a checkpoint without changing model shape, optimizer semantics,
    /// source corpus identity, or deterministic batch sampler.
    public static func resume(model directory: URL, corpus: VivoCellResponseCorpusReader,
                              plan: VivoCellResponseResumePlan, implementation: VivoFingerprint,
                              to destination: URL) throws -> VivoCellResponseModelReceipt {
        return try withGPUExecution {
            try plan.validate()
            try VivoCellResponseArtifactIO.requireNew(destination)
            let loaded = try VivoCellResponseArtifactIO.load(directory, implementation: implementation)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard loaded.state.corpus == corpusFingerprint else {
                throw VivoCellResponseLearningError.incompatible("resume corpus differs from checkpoint")
            }
            let trainedTargetIndices = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            try loaded.trainingPlan.validate(for: corpus)
            let metrics = try runSteps(model: loaded.model, corpus: corpus, plan: loaded.trainingPlan,
                                       startingStep: loaded.state.step, additionalSteps: plan.additionalSteps,
                                       trainedTargetIndices: trainedTargetIndices)
            let state = VivoCellResponseModelState(
                schemaVersion: 4, architecture: loaded.state.architecture, featureAxis: loaded.state.featureAxis,
                targetIDs: loaded.state.targetIDs, trainedTargetBindings: loaded.state.trainedTargetBindings,
                descriptorCount: loaded.state.descriptorCount,
                corpus: loaded.state.corpus, trainingPlan: loaded.state.trainingPlan,
                step: metrics.step, seed: loaded.state.seed, samplerVersion: loaded.state.samplerVersion,
                targetResponsePriorVersion: loaded.state.targetResponsePriorVersion,
                optimizerVersion: loaded.state.optimizerVersion,
                learningRate: loaded.state.learningRate, weightDecay: loaded.state.weightDecay, latestMetrics: metrics)
            return try VivoCellResponseArtifactIO.publish(model: loaded.model, state: state,
                                                           trainingPlanBytes: loaded.trainingPlanBytes,
                                                           implementation: implementation,
                                                           parentCheckpoint: loaded.receipt.checkpoint, to: destination)
        }
    }

    /// Verify receipt, checkpoint framing, safetensors metadata, and every
    /// model parameter shape without opening a source count store.
    public static func verifyModel(_ directory: URL, implementation: VivoFingerprint) throws -> VivoCellResponseModelReceipt {
        try VivoCellResponseArtifactIO.load(directory, implementation: implementation).receipt
    }

    /// Score measured perturbed rows from a held partition. This is a software
    /// evaluation result; it is not a biological validation claim.
    public static func evaluate(model directory: URL, corpus: VivoCellResponseCorpusReader,
                                partition: VivoCellResponsePartition, maximumExamples: Int,
                                implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        return try withGPUExecution {
            let loaded = try VivoCellResponseArtifactIO.load(directory, implementation: implementation)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard loaded.state.corpus == corpusFingerprint else {
                throw VivoCellResponseLearningError.incompatible("evaluation corpus differs from checkpoint")
            }
            let trainedTargetIndices = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            let examples = try corpus.examples(in: partition)
            let result = try evaluate(model: loaded.model, corpus: corpus, examples: examples,
                                      architecture: loaded.state.architecture, trainedTargetIndices: trainedTargetIndices,
                                      maximumExamples: maximumExamples, seed: loaded.state.seed,
                                      salt: 0x4d8b9173)
            let evaluation = VivoCellResponseEvaluation(
                schemaVersion: 2, format: evaluationFormat, partition: partition,
                maximumExamples: maximumExamples, samplerVersion: samplerVersion, seed: loaded.state.seed,
                examples: result.4.count, observedFeatures: result.3,
                negativeLogLikelihood: result.0, rmse: result.1, matchedControlRMSE: result.2,
                selection: result.4)
            try validateEvaluation(evaluation)
            return evaluation
        }
    }

    /// Publish an immutable evaluation receipt after scoring a verified model
    /// and corpus. The stored selection makes later review and replay exact.
    public static func evaluate(model directory: URL, corpus: VivoCellResponseCorpusReader,
                                partition: VivoCellResponsePartition, maximumExamples: Int,
                                implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseEvaluationReceipt {
        let evaluation = try evaluate(model: directory, corpus: corpus, partition: partition,
                                      maximumExamples: maximumExamples, implementation: implementation)
        let modelReceipt = try VivoCellResponseArtifactIO.load(directory, implementation: implementation).receipt
        return try publishEvaluation(evaluation, modelReceipt: modelReceipt, corpus: corpus,
                                     implementation: implementation, to: destination)
    }

    /// Predict a distribution from declared control samples. The output is a
    /// compact response vector artifact, never an H5AD count matrix.
    public static func predict(model directory: URL, corpus: VivoCellResponseCorpusReader,
                               plan: VivoCellResponsePredictionPlan, implementation: VivoFingerprint,
                               to destination: URL) throws -> VivoCellResponsePredictionReceipt {
        return try withGPUExecution {
            let loaded = try VivoCellResponseArtifactIO.load(directory, implementation: implementation)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard loaded.state.corpus == corpusFingerprint else {
                throw VivoCellResponseLearningError.incompatible("prediction corpus differs from checkpoint")
            }
            _ = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            let batch = try inferenceBatch(corpus: corpus, state: loaded.state, plan: plan)
            let output = loaded.model.outputs(context: batch.context, targetIDs: batch.targetIDs,
                                              knownTargetMask: batch.knownTargetMask, descriptors: batch.descriptors)
            let means = output.mean.asArray(Float.self), variances = output.variance.asArray(Float.self),
                baselines = output.baseline.asArray(Float.self), deltas = output.delta.asArray(Float.self)
            guard means.count == corpus.featureCount, variances.count == corpus.featureCount,
                  baselines.count == corpus.featureCount, deltas.count == corpus.featureCount,
                  means.allSatisfy(\.isFinite), baselines.allSatisfy(\.isFinite), deltas.allSatisfy(\.isFinite),
                  variances.allSatisfy({ $0.isFinite && $0 > 0 }) else {
                throw VivoCellResponseLearningError.invalid("prediction output")
            }
            let prediction = VivoCellResponsePrediction(schemaVersion: 1, format: predictionFormat, id: plan.id,
                                                        featureAxis: corpus.plan.featureAxis, featureMask: corpus.featureMask,
                                                        contextBaselineLogCPM: baselines, meanDeltaLogCPM: deltas,
                                                        contextSourceRows: batch.contextRows,
                                                        meanLogCPM: means, varianceLogCPM: variances)
            return try publishPrediction(prediction, query: plan, modelReceipt: loaded.receipt, corpus: corpus,
                                         implementation: implementation, to: destination)
        }
    }

    private static func loadVerifiedPrediction(_ directory: URL,
                                               implementation: VivoFingerprint) throws -> (VivoCellResponsePredictionReceipt, VivoCellResponsePredictionPlan, VivoCellResponsePrediction) {
        let root: VivoRootedDocumentStore
        do {
            root = try .init(rootURL: directory)
        } catch {
            throw VivoCellResponseLearningError.invalid("prediction artifact directory")
        }
        let receiptBytes = try VivoCellResponseArtifactIO.read(root, "receipt.json", maximumBytes: 1_048_576)
        let queryBytes = try VivoCellResponseArtifactIO.read(root, "query.json", maximumBytes: 4_194_304)
        let predictionBytes = try VivoCellResponseArtifactIO.read(root, "prediction.json", maximumBytes: 268_435_456)
        let receipt = try VivoCanonicalJSON.decode(VivoCellResponsePredictionReceipt.self, from: receiptBytes)
        let query = try VivoCanonicalJSON.decode(VivoCellResponsePredictionPlan.self, from: queryBytes)
        let prediction = try VivoCanonicalJSON.decode(VivoCellResponsePrediction.self, from: predictionBytes)
        guard receipt.schemaVersion == 1, receipt.format == predictionFormat, receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == receiptBytes,
              try VivoCanonicalJSON.fingerprint(queryBytes) == receipt.query,
              try VivoCanonicalJSON.fingerprint(predictionBytes) == receipt.prediction,
              try VivoCanonicalJSON.encode(query) == queryBytes,
              try VivoCanonicalJSON.encode(prediction) == predictionBytes,
              prediction.id == query.id else {
            throw VivoCellResponseLearningError.invalid("prediction receipt or output")
        }
        try query.validate()
        try validatePrediction(prediction)
        return (receipt, query, prediction)
    }

    /// Verify a prediction artifact's self-consistency, including its own
    /// query and output fingerprints. This does not establish a model/corpus
    /// provenance chain; use the bound overload when those artifacts exist.
    public static func verifyPrediction(_ directory: URL, implementation: VivoFingerprint) throws -> VivoCellResponsePrediction {
        try loadVerifiedPrediction(directory, implementation: implementation).2
    }

    /// Verify a prediction against the exact model receipt and source-bound
    /// corpus used to create it. This is the production verification route.
    public static func verifyPrediction(_ directory: URL, model modelDirectory: URL,
                                        corpus: VivoCellResponseCorpusReader,
                                        implementation: VivoFingerprint) throws -> VivoCellResponsePrediction {
        return try withGPUExecution {
            let artifact = try loadVerifiedPrediction(directory, implementation: implementation)
            let loaded = try VivoCellResponseArtifactIO.load(modelDirectory, implementation: implementation)
            let modelFingerprint = try VivoCellResponseArtifactIO.modelFingerprint(loaded.receipt)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard artifact.0.model == modelFingerprint,
                  artifact.0.corpus == corpusFingerprint,
                  artifact.2.featureAxis == corpus.plan.featureAxis,
                  artifact.2.featureMask == corpus.featureMask else {
                throw VivoCellResponseLearningError.incompatible("prediction model or corpus provenance differs")
            }
            // Receipt hashes prove stored bytes. Rebuild the exact input and
            // output under the supplied corpus/model as well, so neither an
            // impossible query nor substituted response values can pass.
            _ = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            let batch = try inferenceBatch(corpus: corpus, state: loaded.state, plan: artifact.1)
            let output = loaded.model.outputs(context: batch.context, targetIDs: batch.targetIDs,
                                              knownTargetMask: batch.knownTargetMask, descriptors: batch.descriptors)
            guard artifact.2.contextSourceRows == batch.contextRows,
                  approximatelyEqual(artifact.2.contextBaselineLogCPM, output.baseline.asArray(Float.self)),
                  approximatelyEqual(artifact.2.meanDeltaLogCPM, output.delta.asArray(Float.self)),
                  approximatelyEqual(artifact.2.meanLogCPM, output.mean.asArray(Float.self)),
                  approximatelyEqual(artifact.2.varianceLogCPM, output.variance.asArray(Float.self)) else {
                throw VivoCellResponseLearningError.incompatible("prediction control context or model output differs")
            }
            return artifact.2
        }
    }

    private static func loadVerifiedEvaluation(_ directory: URL,
                                               implementation: VivoFingerprint) throws -> (VivoCellResponseEvaluationReceipt, VivoCellResponseEvaluation) {
        let root: VivoRootedDocumentStore
        do {
            root = try .init(rootURL: directory)
        } catch {
            throw VivoCellResponseLearningError.invalid("evaluation artifact directory")
        }
        let receiptBytes = try VivoCellResponseArtifactIO.read(root, "receipt.json", maximumBytes: 1_048_576)
        let evaluationBytes = try VivoCellResponseArtifactIO.read(root, "evaluation.json", maximumBytes: 134_217_728)
        let receipt = try VivoCanonicalJSON.decode(VivoCellResponseEvaluationReceipt.self, from: receiptBytes)
        let evaluation = try VivoCanonicalJSON.decode(VivoCellResponseEvaluation.self, from: evaluationBytes)
        guard receipt.schemaVersion == 1, receipt.format == evaluationFormat, receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == receiptBytes,
              try VivoCanonicalJSON.encode(evaluation) == evaluationBytes,
              try VivoCanonicalJSON.fingerprint(evaluationBytes) == receipt.evaluation else {
            throw VivoCellResponseLearningError.invalid("evaluation receipt")
        }
        try validateEvaluation(evaluation)
        return (receipt, evaluation)
    }

    /// Verify an evaluation artifact's internal canonical documents. Bind it to
    /// a model and corpus with the overload below before interpreting metrics.
    public static func verifyEvaluation(_ directory: URL, implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        try loadVerifiedEvaluation(directory, implementation: implementation).1
    }

    /// Recompute a published evaluation under the exact model and corpus. This
    /// verifies both its frozen target/control selection and reported metrics.
    public static func verifyEvaluation(_ directory: URL, model modelDirectory: URL,
                                        corpus: VivoCellResponseCorpusReader,
                                        implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        return try withGPUExecution {
            let artifact = try loadVerifiedEvaluation(directory, implementation: implementation)
            let loaded = try VivoCellResponseArtifactIO.load(modelDirectory, implementation: implementation)
            let modelFingerprint = try VivoCellResponseArtifactIO.modelFingerprint(loaded.receipt)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard artifact.0.model == modelFingerprint, artifact.0.corpus == corpusFingerprint,
                  artifact.1.seed == loaded.state.seed else {
                throw VivoCellResponseLearningError.incompatible("evaluation model or corpus provenance differs")
            }
            let trainedTargetIndices = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            let examples = try corpus.examples(in: artifact.1.partition)
            let result = try evaluate(model: loaded.model, corpus: corpus, examples: examples,
                                      architecture: loaded.state.architecture, trainedTargetIndices: trainedTargetIndices,
                                      maximumExamples: artifact.1.maximumExamples, seed: artifact.1.seed,
                                      salt: 0x4d8b9173)
            guard artifact.1.examples == result.4.count, artifact.1.observedFeatures == result.3,
                  artifact.1.selection == result.4,
                  abs(artifact.1.negativeLogLikelihood - result.0) <= 1e-8 * max(1, abs(result.0)),
                  abs(artifact.1.rmse - result.1) <= 1e-8 * max(1, abs(result.1)),
                  abs(artifact.1.matchedControlRMSE - result.2) <= 1e-8 * max(1, abs(result.2)) else {
                throw VivoCellResponseLearningError.incompatible("evaluation selection or metrics differ")
            }
            return artifact.1
        }
    }

    /// Replay an evaluation under its bound model and corpus, then require a
    /// strict improvement over the recorded exact matched-control baseline.
    /// A failed qualification never changes the raw evaluation artifact.
    public static func qualifyEvaluation(_ directory: URL, model modelDirectory: URL,
                                         corpus: VivoCellResponseCorpusReader,
                                         implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        let evaluation = try verifyEvaluation(directory, model: modelDirectory, corpus: corpus,
                                              implementation: implementation)
        try requireMatchedControlImprovement(evaluation)
        return evaluation
    }
}
